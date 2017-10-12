//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// terraces_swath_driver.cpp
// Extract information about terraces using a shapefile of the main stem channel.
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
//
// Fiona Clubb
// University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <ctime>
#include "../LSDRaster.hpp"
#include "../LSDSwathProfile.hpp"
#include "../LSDShapeTools.hpp"
#include "../LSDJunctionNetwork.hpp"
#include "../LSDTerrace.hpp"
#include "../LSDParameterParser.hpp"
#include "../LSDSpatialCSVReader.hpp"

int main (int nNumberofArgs,char *argv[])
{
	//start the clock
	clock_t begin = clock();

	//Test for correct input arguments
  if (nNumberofArgs!=3)
  {
    cout << "=========================================================" << endl;
    cout << "|| Welcome to the terrace swath tool!  	              ||" << endl;
    cout << "|| This program takes in a baseline shapefile and gets ||" << endl;
		cout << "|| a swath profile along a channel from it.            ||" << endl;
		cout << "|| It then extracts the terraces along the channel     ||" << endl;
		cout << "|| using slope and relief thresholds.								  ||" << endl;
    cout << "|| This program was developed by                       ||" << endl;
    cout << "|| Fiona J. Clubb                                      ||" << endl;
    cout << "|| at the University of Edinburgh.                     ||" << endl;
		cout << "|| Was then further developed at the University of     ||" << endl;
		cout << "|| Minnesota, AMERICAN STYLE.                          ||" << endl;
    cout << "=========================================================" << endl;
    cout << "This program requires two inputs: " << endl;
    cout << "* First the path to the parameter file." << endl;
    cout << "* Second the name of the param file (see below)." << endl;
    cout << "---------------------------------------------------------" << endl;
    cout << "Then the command line argument will be, for example: " << endl;
    cout << "In linux:" << endl;
    cout << "./get_terraces_from_shapefile.out /LSDTopoTools/Topographic_projects/Test_data/ LSDTT_terraces.param" << endl;
		cout << "For more information please see the documentation: " << endl;
		cout << "http://lsdtopotools.github.io/LSDTT_book/#_terraces" << endl;
		cout << "Although the documentation is not up to date for this version, SORRY." << endl;
    cout << "=========================================================" << endl;
    exit(EXIT_SUCCESS);
  }
  string path_name = argv[1];
  string f_name = argv[2];

	// maps for setting default parameters
	map<string,int> int_default_map;
	map<string,float> float_default_map;
	map<string,bool> bool_default_map;
	map<string,string> string_default_map;

	// set default int parameters
	int_default_map["Threshold_SO"] = 3;
	int_default_map["Relief lower percentile"] = 25;
	int_default_map["Relief upper percentile"] = 75;
	int_default_map["Slope lower percentile"] = 25;
	int_default_map["Slope upper percentile"] = 75;
	int_default_map["Min patch size"] = 1000;
	int_default_map["search_radius"] = 10;
	int_default_map["NormaliseToBaseline"] = 1;
	int_default_map["Min terrace height"] = 2;
	int_default_map["Chan area threshold"] = 1000;

	// set default float parameters
	float_default_map["surface_fitting_window_radius"] = 6;
	float_default_map["Min slope filling"] = 0.0001;
	float_default_map["QQ threshold"] = 0.005;
	float_default_map["HalfWidth"] = 500;

	// set default bool parameters
	bool_default_map["Filter topography"] = true;

	// set default string parameters
	string_default_map["input_shapefile"] = "NULL";

	// Use the parameter parser to get the maps of the parameters required for the
	// analysis
	// load parameter parser object
	LSDParameterParser LSDPP(path_name,f_name);
	LSDPP.force_bil_extension();

	LSDPP.parse_all_parameters(float_default_map, int_default_map, bool_default_map,string_default_map);
	map<string,float> this_float_map = LSDPP.get_float_parameters();
	map<string,int> this_int_map = LSDPP.get_int_parameters();
	map<string,bool> this_bool_map = LSDPP.get_bool_parameters();
	map<string,string> this_string_map = LSDPP.get_string_parameters();

	// Now print the parameters for bug checking
	LSDPP.print_parameters();

	// location of the files
	string DATA_DIR =  LSDPP.get_read_path();
	string DEM_ID =  LSDPP.get_read_fname();
	string OUT_DIR = LSDPP.get_write_path();
	string OUT_ID = LSDPP.get_write_fname();
	string DEM_extension =  LSDPP.get_dem_read_extension();
	vector<string> boundary_conditions = LSDPP.get_boundary_conditions();
	string CHeads_file = LSDPP.get_CHeads_file();

	if (this_string_map["coords_csv_file"] == "NULL")
	{
		cout << "FATAL ERROR: I can't find your coordinates file. Check your spelling!! \n The parameter key needs to be 'coords_csv_file'" << endl;
		exit(EXIT_SUCCESS);
	}

  cout << "starting the test run... here we go!" << endl;

	LSDRaster RasterTemplate;

	if(this_bool_map["Filter topography"])
	{
		 // load the DEM
		 cout << "Loading the DEM..." << endl;
		 LSDRaster load_DEM((DATA_DIR+DEM_ID), DEM_extension);
		 RasterTemplate = load_DEM;

		 // filter using Perona Malik
		 int timesteps = 50;
		 float percentile_for_lambda = 90;
		 float dt = 0.1;
		 RasterTemplate = RasterTemplate.PeronaMalikFilter(timesteps, percentile_for_lambda, dt);

		 // fill
		 RasterTemplate = RasterTemplate.fill(this_float_map["Min slope filling"]);
		 string fill_name = "_filtered";
		 RasterTemplate.write_raster((DATA_DIR+DEM_ID+fill_name), DEM_extension);
	}
	else
	{
		//don't do the filtering, just load the filled DEM
		LSDRaster load_DEM((DATA_DIR+DEM_ID+"_filtered"), DEM_extension);
		RasterTemplate = load_DEM;
	}

	cout << "\t Flow routing..." << endl;
	// get a flow info object
	LSDFlowInfo FlowInfo(boundary_conditions, RasterTemplate);
	// calcualte the distance from outlet
	LSDRaster DistanceFromOutlet = FlowInfo.distance_from_outlet();

	// some error checking
	vector<int> sources;
	if (CHeads_file == "NULL")
	{
		cout << "I can't find your channel heads file so I'm going to use an area threshold to extract the sources" << endl;
		LSDIndexRaster ContributingPixels = FlowInfo.write_NContributingNodes_to_LSDIndexRaster();
		sources = FlowInfo.get_sources_index_threshold(ContributingPixels, this_int_map["Chan_area_threshold"]);
	}
	else
	{
		cout << "\t Loading the sources" << endl;
		// load the sources
		vector<int> sources = FlowInfo.Ingest_Channel_Heads((DATA_DIR+CHeads_file), "csv", 2);
		cout << "\t Got sources!" << endl;
	}

	// now get the junction network
	LSDJunctionNetwork ChanNetwork(sources, FlowInfo);
  cout << "\t Got the channel network" << endl;

	cout << "\t loading baseline points" << endl;
	PointData BaselinePoints = LoadShapefile(path_name+this_string_map["input_shapefile"].c_str());

  cout << "\t Creating swath template" << endl;
  LSDSwath TestSwath(BaselinePoints, RasterTemplate, this_float_map["HalfWidth"]);

	cout << "\n\t Getting raster from swath" << endl;
	LSDRaster SwathRaster = TestSwath.get_raster_from_swath_profile(RasterTemplate, this_int_map["NormaliseToBaseline"]);
	string swath_ext = "_swath_raster";
	SwathRaster.write_raster((DATA_DIR+DEM_ID+swath_ext), DEM_extension);

	// get the elevation raster from the swath
	//LSDRaster ElevationRaster = TestSwath.get_raster_from_swath_profile(RasterTemplate, 0);

  // get the slope
	cout << "\t Getting the slope" << endl;
  vector<LSDRaster> surface_fitting;
  LSDRaster Slope;
  vector<int> raster_selection(8, 0);
  raster_selection[1] = 1;             // this means you want the slope
  surface_fitting = RasterTemplate.calculate_polyfit_surface_metrics(this_float_map["surface_fitting_window_radius"], raster_selection);
  Slope = surface_fitting[1];

	float mask_threshold = 1.0;
	bool below = 0;
	// remove any stupid slope values
	LSDRaster Slope_new = Slope.mask_to_nodata_using_threshold(mask_threshold, below);

	// get the channel relief and slope threshold using quantile-quantile plots
	cout << "Getting channel relief threshold from QQ plots" << endl;
	string qq_fname = DATA_DIR+DEM_ID+"_qq_relief.txt";
	float relief_threshold_from_qq = SwathRaster.get_threshold_for_floodplain_QQ(qq_fname, this_float_map["QQ threshold"], this_int_map["Relief lower percentile"], this_int_map["Relief upper percentile"]);

	cout << "Getting slope threshold from QQ plots" << endl;
	string qq_slope = DATA_DIR+DEM_ID+"_qq_slope.txt";
	float slope_threshold_from_qq = Slope_new.get_threshold_for_floodplain_QQ(qq_slope, this_float_map["QQ threshold"], this_int_map["Slope lower percentile"], this_int_map["Slope upper percentile"]);

	cout << "Relief threshold: " << relief_threshold_from_qq << " Slope threshold: " << slope_threshold_from_qq << endl;

	// get the terrace pixels
	LSDTerrace Terraces(SwathRaster, Slope_new, ChanNetwork, FlowInfo, relief_threshold_from_qq, slope_threshold_from_qq, this_int_map["Min patch size"], this_int_map["Threshold_SO"], this_int_map["Min terrace height"]);
	LSDIndexRaster ConnectedComponents = Terraces.print_ConnectedComponents_to_Raster();
	string CC_ext = "_terrace_IDs";
	ConnectedComponents.write_raster((DATA_DIR+DEM_ID+CC_ext), DEM_extension);

	// cout << "\t Testing connected components" << endl;
	// vector <vector <float> > CC_vector = TestSwath.get_connected_components_along_swath(ConnectedComponents, RasterTemplate, this_int_map["NormaliseToBaseline"]);
	//
	// // push back results to file for plotting
	// ofstream output_file_CC;
	// string output_fname = "_terrace_swath_plots.txt";
	// output_file_CC.open((DATA_DIR+DEM_ID+output_fname).c_str());
	// for (int i = 0; i < int(CC_vector[0].size()); ++i)
	// {
	// 	output_file_CC << CC_vector[0][i] << " " << CC_vector[1][i] << " " << CC_vector[2][i] << endl;
	// }
	// output_file_CC.close();

	// write raster of terrace elevations
	LSDRaster ChannelRelief = Terraces.get_Terraces_RasterValues(SwathRaster);
	string relief_ext = "_terrace_relief_final";
	ChannelRelief.write_raster((DATA_DIR+DEM_ID+relief_ext), DEM_extension);

	// print the terrace information to a csv
	string csv_fname = "_terrace_info.csv";
	string full_csv_name = DATA_DIR+DEM_ID+csv_fname;
	cout << "The full csv filename is: " << full_csv_name << endl;
	Terraces.print_TerraceInfo_to_csv(full_csv_name, RasterTemplate, ChannelRelief, FlowInfo, TestSwath);

	// Done, check how long it took
	clock_t end = clock();
	float elapsed_secs = float(end - begin) / CLOCKS_PER_SEC;
	cout << "DONE, Time taken (secs): " << elapsed_secs << endl;
}
