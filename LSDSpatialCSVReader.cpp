//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// LSDSpatialCSVReader.hpp
// Land Surface Dynamics SpatialCSVReader
//
// An object within the University
//  of Edinburgh Land Surface Dynamics group topographic toolbox
//  for reading csv data. The data needs to have latitude and longitude
//  in WGS84 coordinates.
//
// Developed by:
//  Simon M. Mudd
//  Martin D. Hurst
//  David T. Milodowski
//  Stuart W.D. Grieve
//  Declan A. Valters
//  Fiona Clubb
//
// Copyright (C) 2017 Simon M. Mudd 2017
//
// Developer can be contacted by simon.m.mudd _at_ ed.ac.uk
//
//    Simon Mudd
//    University of Edinburgh
//    School of GeoSciences
//    Drummond Street
//    Edinburgh, EH8 9XP
//    Scotland
//    United Kingdom
//
// This program is free software;
// you can redistribute it and/or modify it under the terms of the
// GNU General Public License as published by the Free Software Foundation;
// either version 2 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY;
// without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the
// GNU General Public License along with this program;
// if not, write to:
// Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor,
// Boston, MA 02110-1301
// USA
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
#include <fstream>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <ctype.h>
#include <sstream>
#include <algorithm>
#include <vector>
#include "LSDStatsTools.hpp"
#include "LSDShapeTools.hpp"
#include "LSDCosmoData.hpp"
#include "LSDRaster.hpp"
#include "LSDFlowInfo.hpp"
#include "LSDJunctionNetwork.hpp"
#include "LSDBasin.hpp"
#include "LSDSpatialCSVReader.hpp"
#include "LSDRasterInfo.hpp"
#include "TNT/tnt.h"
using namespace std;
using namespace TNT;

#ifndef LSDSpatialCSVReader_CPP
#define LSDSpatialCSVReader_CPP



// empty create function
void LSDSpatialCSVReader::create()
{
  cout << "Size data map: " << data_map.size() << endl;
  map<string, vector<string> > empty_map;
  cout << "Size empty map: " << empty_map.size() << endl;
  data_map = empty_map;
}

//==============================================================================
// Basic create function
//==============================================================================
void LSDSpatialCSVReader::create(string csv_fname)
{
  cout << "I am creating something for you" << endl;
  NRows = -9999;
  NCols = -9999;
  XMinimum = -9999;
  YMinimum = -9999;
  DataResolution = -9999;
  NoDataValue = -9999;

  ///A map of strings for holding georeferencing information
  map<string,string> EmptyString;
  GeoReferencingStrings = EmptyString;

  //cout << "Size data map: " << data_map.size() << endl;
  //cout << "Size test map: " << test_map.size() << endl;

  //map<string, vector<string> > empty_map;
  //cout << "Size empty map: " << empty_map.size() << endl;
  //data_map = empty_map;


  load_csv_data(csv_fname);

}



//==============================================================================
// Basic create function
//==============================================================================
void LSDSpatialCSVReader::create(LSDRasterInfo& ThisRasterInfo, string csv_fname)
{
  cout << "I am creating something for you" << endl;
  NRows = ThisRasterInfo.get_NRows();
  NCols = ThisRasterInfo.get_NCols();
  XMinimum = ThisRasterInfo.get_XMinimum();
  YMinimum = ThisRasterInfo.get_YMinimum();
  DataResolution = ThisRasterInfo.get_DataResolution();
  NoDataValue = ThisRasterInfo.get_NoDataValue();
  GeoReferencingStrings = ThisRasterInfo.get_GeoReferencingStrings();

  //cout << "Size data map: " << data_map.size() << endl;
  //cout << "Size test map: " << test_map.size() << endl;

  //map<string, vector<string> > empty_map;
  //cout << "Size empty map: " << empty_map.size() << endl;
  //data_map = empty_map;


  load_csv_data(csv_fname);

}

//==============================================================================
// Basic create function
//==============================================================================
void LSDSpatialCSVReader::create(LSDRaster& ThisRaster, string csv_fname)
{
  cout << "I am creating something for you2" << endl;
  NRows = ThisRaster.get_NRows();
  NCols = ThisRaster.get_NCols();
  XMinimum = ThisRaster.get_XMinimum();
  YMinimum = ThisRaster.get_YMinimum();
  DataResolution = ThisRaster.get_DataResolution();
  NoDataValue = ThisRaster.get_NoDataValue();
  GeoReferencingStrings = ThisRaster.get_GeoReferencingStrings();

  load_csv_data(csv_fname);
}


//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This loads a csv file
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSpatialCSVReader::load_csv_data(string filename)
{
  // make sure the filename works
  ifstream ifs(filename.c_str());
  if( ifs.fail() )
  {
    cout << "\nFATAL ERROR: Trying to load csv cosmo data file, but the file" << filename
         << "doesn't exist; LINE 245 LSDCosmoData" << endl;
    exit(EXIT_FAILURE);
  }
  else
  {
    cout << "I have opened the csv file." << endl;
  }

  // Initiate the data map
  map<string, int > temp_vec_vec_key;
  vector< vector<string> > temp_vec_vec;
  map<string, vector<string> > temp_data_map;

  cout << "Data map size is: " << data_map.size() << endl;
  cout << "longitude size is: " << longitude.size() << endl;
  data_map = temp_data_map;

  // initiate the string to hold the file
  string line_from_file;
  vector<string> empty_string_vec;
  vector<string> this_string_vec;
  string temp_string;

  vector<double> temp_longitude;
  vector<double> temp_latitude;

  // get the headers from the first line
  getline(ifs, line_from_file);

  // reset the string vec
  this_string_vec = empty_string_vec;

  // create a stringstream
  stringstream ss(line_from_file);
  ss.precision(9);

  while( ss.good() )
  {
    string substr;
    getline( ss, substr, ',' );

    // remove the spaces
    substr.erase(remove_if(substr.begin(), substr.end(), ::isspace), substr.end());

    // remove control characters
    substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());

    // add the string to the string vec
    this_string_vec.push_back( substr );
  }
  // now check the data map
  int n_headers = int(this_string_vec.size());
  vector<string> header_vector = this_string_vec;
  int latitude_index = -9999;
  int longitude_index = -9999;
  for (int i = 0; i<n_headers; i++)
  {
    //cout << "This header is: " << this_string_vec[i] << endl;
    if (this_string_vec[i]== "latitude")
    {
      latitude_index = i;
      //cout << "The latitude index is: " << latitude_index << endl;

    }
    else if (this_string_vec[i] == "longitude")
    {
      longitude_index = i;
      //cout << "The longitude index is: " << longitude_index << endl;
    }
    else
    {
      temp_data_map[header_vector[i]] = empty_string_vec;
    }
  }


  // now loop through the rest of the lines, getting the data.
  while( getline(ifs, line_from_file))
  {
    //cout << "Getting line, it is: " << line_from_file << endl;
    // reset the string vec
    this_string_vec = empty_string_vec;

    // create a stringstream
    stringstream ss(line_from_file);

    while( ss.good() )
    {
      string substr;
      getline( ss, substr, ',' );

      // remove the spaces
      substr.erase(remove_if(substr.begin(), substr.end(), ::isspace), substr.end());

      // remove control characters
      substr.erase(remove_if(substr.begin(), substr.end(), ::iscntrl), substr.end());

      // add the string to the string vec
      this_string_vec.push_back( substr );
    }

    //cout << "Yoyoma! size of the string vec: " <<  this_string_vec.size() << endl;
    if ( int(this_string_vec.size()) <= 0)
    {
      cout << "Hey there, I am trying to load your csv data but you seem not to have" << endl;
      cout << "enough columns in your file. I am ignoring a line" << endl;
    }
    else
    {
      int n_cols = int(this_string_vec.size());
      //cout << "N cols is: " << n_cols << endl;
      for (int i = 0; i<n_cols; i++)
      {
        if (i == latitude_index)
        {
          temp_latitude.push_back( atof(this_string_vec[i].c_str() ) );
        }
        else if (i == longitude_index)
        {
          temp_longitude.push_back( atof(this_string_vec[i].c_str() ) );

        }
        else
        {
          temp_data_map[header_vector[i]].push_back(this_string_vec[i]);
        }

      }
      //cout << "Done with this line." << endl;
    }

  }



  //cout << "Assigning the vectors." << endl;
  latitude = temp_latitude;
  longitude = temp_longitude;
  data_map = temp_data_map;
  //cout << "Done reading your file." << endl;
}
//==============================================================================



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This returns the string vector of data from a given column name
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<string> LSDSpatialCSVReader::get_data_column(string column_name)
{
  vector<string> data_vector;
  if ( data_map.find(column_name) == data_map.end() )
  {
    // not found
    cout << "I'm afraid the column "<< column_name << " is not in this dataset" << endl;
  }
  else
  {
    data_vector = data_map[column_name];
  }
  return data_vector;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function returns data that is in the raster for snapping
// It focuses on the ID vector, and converts to Easting-Northing coordinates
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSpatialCSVReader::get_data_in_raster_for_snapping(string column_name,
                                    vector<float>& UTMEasting,
                                    vector<float>& UTMNorthing,
                                    vector<string>& data_vector)
{
  vector<string> this_data_vector;
  vector<string> thinned_data_vector;
  vector<float> easting;
  vector<float> northing;
  if ( data_map.find(column_name) == data_map.end() )
  {
    // not found
    cout << "I am afraid you tried to access data in the csv file that isn't there." << endl;
    cout << "The missing column name is: " << column_name << endl;
    exit(EXIT_SUCCESS);
  }
  else
  {
    this_data_vector = data_map[column_name];

  }

  // make sure the vector of booleans stating if the point is in the DEM exists
  check_if_points_are_in_raster();

  // convert to UTM
  vector<float> UTME;
  vector<float> UTMN;
  get_x_and_y_from_latlong(UTME,UTMN);

  // now loop through the samples
  int N_samples = int(longitude.size());
  for (int i = 0; i<N_samples; i++)
  {
    if (is_point_in_raster[i])
    {
      easting.push_back(UTME[i]);
      northing.push_back(UTMN[i]);
      thinned_data_vector.push_back(this_data_vector[i]);

      cout << "I am pushing E: " << UTME[i] << " N: " << UTMN[i] << " and data: " << this_data_vector[i] << endl;

    }
  }

  UTMEasting = easting;
  UTMNorthing = northing;
  data_vector = thinned_data_vector;
}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-



//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// Converts a data column to a float vector
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
vector<float> LSDSpatialCSVReader::data_column_to_float(string column_name)
{
  vector<string> string_vec = get_data_column(column_name);
  vector<float> float_vec;
  int N_data_elements = string_vec.size();
  for(int i = 0; i<N_data_elements; i++)
  {
    float_vec.push_back( atof(string_vec[i].c_str()));
  }
  return float_vec;
}

// Converts a data column to a float vector
vector<int> LSDSpatialCSVReader::data_column_to_int(string column_name)
{
  vector<string> string_vec = get_data_column(column_name);
  vector<int> int_vec;
  int N_data_elements = string_vec.size();
  if (N_data_elements == 0)
  {
    cout << "Couldn't read in the data column. Check the column name!" << endl;
  }
  for(int i = 0; i<N_data_elements; i++)
  {
    int_vec.push_back( atoi(string_vec[i].c_str()));
  }
  return int_vec;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This function gets the UTM zone
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSpatialCSVReader::get_UTM_information(int& UTM_zone, bool& is_North)
{

  // set up strings and iterators
  map<string,string>::iterator iter;

  //check to see if there is already a map info string
  string mi_key = "ENVI_map_info";
  iter = GeoReferencingStrings.find(mi_key);

  if (iter != GeoReferencingStrings.end() )
  {
    string info_str = GeoReferencingStrings[mi_key] ;

    // now parse the string
    vector<string> mapinfo_strings;
    istringstream iss(info_str);
    while( iss.good() )
    {
      string substr;
      getline( iss, substr, ',' );
      mapinfo_strings.push_back( substr );
    }
    UTM_zone = atoi(mapinfo_strings[7].c_str());
    //cout << "Line 1041, UTM zone: " << UTM_zone << endl;
    //cout << "LINE 1042 LSDRaster, N or S: " << mapinfo_strings[7] << endl;

    // find if the zone is in the north
    string n_str = "n";
    string N_str = "N";
    is_North = false;
    size_t found = mapinfo_strings[8].find(N_str);
    if (found!=std::string::npos)
    {
      is_North = true;
    }
    found = mapinfo_strings[8].find(n_str);
    if (found!=std::string::npos)
    {
      is_North = true;
    }
    //cout << "is_North is: " << is_North << endl;

  }
  else
  {
    UTM_zone = NoDataValue;
    is_North = false;
  }

}
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// This sets some coordinate system strings for UTM
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSpatialCSVReader::set_UTM_information(int UTM_zone, bool is_North)
{
  string cs_key = "coordinate system string";

  string fpart = "{PROJCS[\"WGS_1984_UTM_Zone_";
  string spart = itoa(UTM_zone);

  string tpart;
  if(is_North)
  {
    tpart = "N";
  }
  else
  {
    tpart = "S";
  }
  string fopart = ",GEOGCS[\"GCS_WGS_1984\",DATUM[\"D_WGS_1984\",SPHEROID[\"WGS_1984\",6378137,298.257223563]],PRIMEM[\"Greenwich\",0],UNIT[\"Degree\",0.017453292519943295]],PROJECTION[\"Transverse_Mercator\"],PARAMETER[\"latitude_of_origin\",0],PARAMETER[\"central_meridian\",";

  int start_cm = -183;
  int cm = start_cm+UTM_zone*6;
  string fipart = itoa(cm);
  string sipart = "],PARAMETER[\"scale_factor\",0.9996],PARAMETER[\"false_easting\",500000],PARAMETER[\"false_northing\",0],UNIT[\"Meter\",1]]}";

  string cs_string = fpart+spart+tpart+fopart+fipart+sipart;
  GeoReferencingStrings[cs_key] = cs_string;
  cout << "The coordinate system string is: " << cs_string << endl;

  string mi_key = "ENVI_map_info";
  string mi_p1 = "{UTM, 1, 1, -9999, -9999, -9999, -9999,";
  string mi_p2 = itoa(UTM_zone);

  string mi_p3;
  if(is_North)
  {
    mi_p3 = ", North,WGS-84}";
  }
  else
  {
    mi_p3 = ", South,WGS-84}";
  }

  string mi = mi_p1+mi_p2+mi_p3;
  GeoReferencingStrings[mi_key] = mi;
  cout << "the map info string is:" << mi << endl;

}


//==============================================================================
// This gets the x and y locations from the latitude and longitude
//==============================================================================
void LSDSpatialCSVReader::get_x_and_y_from_latlong(vector<float>& UTME,vector<float>& UTMN)
{
  // initilise the converter
  LSDCoordinateConverterLLandUTM Converter;

  int N_samples =  int(latitude.size());

  // set up some temporary vectors
  vector<float> this_UTMN(N_samples,0);
  vector<float> this_UTME(N_samples,0);

  double this_Northing;
  double this_Easting;

  int UTM_zone;
  bool is_North;
  get_UTM_information(UTM_zone, is_North);


  // loop throught the samples collecting UTM information
  int eId = 22;             // defines the ellipsiod. This is WGS
  for(int i = 0; i<N_samples; i++)
  {
    //cout << "Converting point " << i << " to UTM." << endl;
    Converter.LLtoUTM_ForceZone(eId, latitude[i], longitude[i],
                      this_Northing, this_Easting, UTM_zone);
    this_UTMN[i] = this_Northing;
    this_UTME[i] = this_Easting;
    //cout << "Easting: " << this_Easting << " and northing: " << this_Northing << endl;
  }

  UTME = this_UTME;
  UTMN = this_UTMN;
}

//==============================================================================
// This gets the latitude and longitude from x and y columns
//==============================================================================
void LSDSpatialCSVReader::get_latlong_from_x_and_y(string X_column_name, string Y_column_name)
{
  // initilise the converter
  LSDCoordinateConverterLLandUTM Converter;

  vector<double> new_lat;
  vector<double> new_long;

  int UTM_zone;
  bool is_North;
  int eId = 22;
  get_UTM_information(UTM_zone,is_North);
  cout << "Getting lat and long from UTM Easting and Northing. " << endl;
  cout << "Zone: " << UTM_zone << " and is north? ";
  if (is_North)
  {
    cout <<  "youbetcha!" << endl;
  }
  else
  {
    cout << " no, it is south." << endl;
  }

  vector<float> X_data = data_column_to_float(X_column_name);
  vector<float> Y_data = data_column_to_float(Y_column_name);

  int N_x = int(X_data.size());
  int N_y = int(Y_data.size());
  if (N_x != N_y)
  {
    cout << "Your X and Y columns don't have the same lengths, something has gone wrong." << endl;
    cout << "I am not updating the latitude and longitude" << endl;
  }
  else
  {
    for (int i = 0; i<N_x; i++)
    {
      double thisX = X_data[i];
      double thisY = Y_data[i];

      double Lat;
      double Long;
      Converter.UTMtoLL(eId, thisY, thisX, UTM_zone, is_North,Lat, Long);
      new_lat.push_back(Lat);
      new_long.push_back(Long);

      //if (i == 0)
      //{
      //  cout << "X: " << thisX << " Y: " << thisY << " Lat: " << Lat << " Long: " << Long << endl;
      //}

    }

    latitude = new_lat;
    longitude = new_long;
  }

}


//==============================================================================
// This checks if points are in raster
//==============================================================================
void LSDSpatialCSVReader::check_if_points_are_in_raster()
{
  vector<float> UTME;   // easting
  vector<float> UTMN;   // northing
  vector<bool> temp_is_point_in_raster;

  // get the easting and northing
  get_x_and_y_from_latlong(UTME,UTMN);

  bool is_in_raster;

  int N_samples = int(latitude.size());
  for(int i = 0; i<N_samples; i++)
  {

    is_in_raster = true;

    // Shift origin to that of dataset
    float X_coordinate_shifted_origin = UTME[i] - XMinimum;
    float Y_coordinate_shifted_origin = UTMN[i] - YMinimum;

    // Get row and column of point
    int col_point = int(X_coordinate_shifted_origin/DataResolution);
    int row_point = (NRows - 1) - int(round(Y_coordinate_shifted_origin/DataResolution));

    if(col_point < 0 || col_point > NCols-1 || row_point < 0 || row_point > NRows -1)
    {
      is_in_raster = false;
    }
    temp_is_point_in_raster.push_back(is_in_raster);
  }

  is_point_in_raster = temp_is_point_in_raster;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Function to get vectors of x and y coordinates, and the node indices of these
// points - this DOES NOT use latitude and longitude, instead it assumes that
// your csv file has columns labelled "X" and "Y". Can be used when you want to
// read in points without converting from latitude and longitude.
// FJC 21/02/17
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSpatialCSVReader::get_nodeindices_from_x_and_y_coords(LSDFlowInfo& FlowInfo, vector<float>& X_coords, vector<float>& Y_coords, vector<int> & NodeIndices)
{
  // get vectors from the columns in the csv file
  string X_coords_name = "X";
  string Y_coords_name = "Y";
  vector<float> X_coords_temp = data_column_to_float(X_coords_name);
  vector<float> Y_coords_temp = data_column_to_float(Y_coords_name);
  vector<int> NodeIndices_temp;

  for (int i = 0; i < int(X_coords_temp.size()); ++i)
  {
    int NodeIndex = FlowInfo.get_node_index_of_coordinate_point(X_coords_temp[i], Y_coords_temp[i]);
    if (NodeIndex != NoDataValue) { NodeIndices_temp.push_back(NodeIndex); }
  }

  //copy to output vectors
  X_coords = X_coords_temp;
  Y_coords = Y_coords_temp;
  NodeIndices = NodeIndices_temp;
}

//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This prints the data column keys to screen
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
void LSDSpatialCSVReader::print_data_map_keys_to_screen()
{
  cout << "These are the keys: " << endl;
  for( map<string, vector<string> >::iterator it = data_map.begin(); it != data_map.end(); ++it)
  {
    cout << "Key is: " <<it->first << "\n";
  }
}

//==============================================================================
// This prints the lat and long to screen
//==============================================================================
void LSDSpatialCSVReader::print_lat_long_to_screen()
{
  int N_data = int(latitude.size());
  cout << "latitude,longitude"<< endl;
  cout.precision(9);
  for (int i = 0; i< N_data; i++)
  {
    cout << latitude[i] << "," << longitude[i] << endl;
  }
}

//==============================================================================
// This prints the lat and long to screen
//==============================================================================
void LSDSpatialCSVReader::print_lat_long_to_screen(bool only_print_in_raster)
{
  if (only_print_in_raster)
  {
    check_if_points_are_in_raster();
  }


  int N_data = int(latitude.size());
  cout << "latitude,longitude"<< endl;
  cout.precision(9);
  for (int i = 0; i< N_data; i++)
  {
    if (only_print_in_raster)
    {
      if (is_point_in_raster[i])
      {
        cout << latitude[i] << "," << longitude[i] << endl;
      }

    }
    else
    {
      cout << latitude[i] << "," << longitude[i] << endl;
    }
  }
}

//==============================================================================
// This prints the UTM coordinates to csv for checking
// FJC 03/03/17
//==============================================================================
void LSDSpatialCSVReader::print_UTM_coords_to_csv(vector<float> UTME, vector<float> UTMN, string csv_outname)
{
  ofstream outfile;
  outfile.open(csv_outname.c_str());

  outfile << "PointNo,x,y" << endl;
  for (int i = 0; i < int(UTME.size()); i++)
  {
    outfile << i+1 << "," << UTME[i] << "," << UTMN[i] << endl;
  }
  outfile.close();
}


#endif
