//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
// LSDSwathProfile.hpp
//------------------------------------------------------------------------------
// This code houses the LSDSwath object, used to make swath profiles
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// This object is written by
// David T. Milodowski, University of Edinburgh
// Simon M. Mudd, University of Edinburgh
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
//
// Version 0.0.1		17/02/2014
// Prerequisite software packages: TNT, PCL and liblas
//
//=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

#include <string>
#include <vector>
#include "TNT/tnt.h"
#include "LSDRaster.hpp"
#include "LSDCloudRaster.hpp"
#include "LSDShapeTools.hpp"
#include "LSDIndexChannel.hpp"
// PCL
#include <pcl/point_types.h>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>
#include <pcl/octree/octree.h>
// liblas
//#include <liblas/liblas.hpp>
using namespace std;
using namespace TNT;

#ifndef LSDSwathProfile_H
#define LSDSwathProfile_H

/// @brief This code houses the LSDSwath object, used to make swath profiles
/// @author DTM
/// @date 17/02/14
class LSDSwath
{
  public:
  LSDSwath()	{ create(); }
  ///@brief create an LSDSwath using a raster as a template.
  ///
  ///@param PointData ProfilePoints -> coordinates of points forming the
  /// baseline of the profile.
  ///@param LSDRaster RasterTemplate -> a raster dataset that is used as a
  /// template for the swath profile i.e. any rasters that you wish to generate
  /// the profile for should have the same characteristics/spatial extent as the
  /// original template.
  ///@param float ProfileHalfWidth
  ///@author DTM
  ///@date 11/04/2014
  ///
  LSDSwath(PointData& ProfilePoints, LSDRaster& RasterTemplate, float HalfWidth) { create(ProfilePoints, RasterTemplate, HalfWidth); }

  ///@brief create an LSDSwath using a raster as a template.
  ///
  ///@param vector<vector<float> >& Y_X_points -> coordinates of points
  /// forming the extremities of the baseline of the profile.
  ///@param LSDRaster RasterTemplate -> a raster dataset that is used as a
  /// template for the swath profile i.e. any rasters that you wish to generate
  /// the profile for should have the same characteristics/spatial extent as the
  /// original template.
  ///@param float ProfileHalfWidth
  ///@author DTM
  ///@date 28/02/2017
  ///
  LSDSwath(vector<float>& Y_X_points, LSDRaster& RasterTemplate, float& HalfWidth, float d_space) { create(Y_X_points, RasterTemplate, HalfWidth, d_space); }

  ///@brief create an LSDSwath from a series of csv files
  ///@param path the path name
  ///@param RasterTemplate a lsd raster template, needs to be the same as
  /// the one originally used to create the swath
  ///@param DEM_prefix name of the DEM_prefix
  ///@Param FlowInfo LSDFlowInfo object
  ///@author FJC
  ///@date 13/11/17
  LSDSwath(string path, LSDRaster& RasterTemplate, string DEM_prefix, LSDFlowInfo& FlowInfo) { create(path, RasterTemplate, DEM_prefix, FlowInfo); }

  void get_transverse_swath_profile(LSDRaster& Raster, vector<float> desired_percentiles, float BinWidth,
       vector<float>& mid_points, vector<float>& mean_profile, vector<float>& sd_profile, vector< vector<float> >& output_percentile_profiles,
       int NormaliseToBaseline);

  void get_longitudinal_swath_profile(LSDRaster& Raster, vector<float> desired_percentiles, float BinWidth,
       vector<float>& mid_points, vector<float>& mean_profile, vector<float>& sd_profile, vector< vector<float> >& output_percentile_profiles,
       int NormaliseToBaseline);

  ///@brief create a raster in the shape of the swath profile
  ///@param Raster LSDRaster of interest
  ///@param NormaliseToBaseline if 0 --> raster values; if 1 --> raster values normalised to baseline value
  ///@return raster in shape of swath profile
  ///@author FJC
  ///@date 16/10/15
  LSDRaster get_raster_from_swath_profile(LSDRaster& Raster, int NormaliseToBaseline);

	///@brief fill in the baseline raster value with the average value of the pixels along the transverse swath profile
	///@param Raster Raster template
	///@return LSDRaster with filled in values along the baseline
	///@author FJC
	///@date 16/01/17
	LSDRaster fill_in_channels_swath(LSDRaster& Raster);

  ///@brief Get information about connected components along the swath profile
  ///@details This function takes in a connected components raster and returns the average
  /// value of another chosen raster and distance along the baseline of each
  /// connected components patch. The user can choose whether to normalise the
  /// second raster to the baseline value.
  /// vector of vectors has the format:
  /// 0 = patch ID
  /// 1 = mean raster value for the patch id
  /// 2 = mean distance along the baseline
  ///@return vector of vector with patch ids, mean values, and distance along baseline
  ///@author FJC
  ///@date 24/01/17
  vector <vector <float> > get_connected_components_along_swath(LSDIndexRaster& ConnectedComponents, LSDRaster& RasterTemplate, int NormaliseToBaseline);

  ///@details This function takes in a raster and returns the mean, min and max values of the raster
  /// at each point along the swath
  /// vector of vectors has the format:
  /// 0 = distance along swath
  /// 1 = mean value along swath
  /// 2 = min value along swath
  /// 3 = max value along swath
  /// if NormaliseToBaseline == 1 then the values will be normalised to the baseline.
  ///@return vector of vectors
  ///@author FJC
  ///@date 15/02/17
  vector <vector <float> > get_RasterValues_along_swath(LSDRaster& RasterTemplate, int NormaliseToBaseline);

  ///@details This function takes in a connected components raster and returns an array
  /// of the distance along the baseline of each pixel in the raster
  ///@param ConnectedComponents connected components raster
  ///@return array with baseline components
  ///@author FJC
  ///@date 28/09/17
  Array2D<float> get_BaselineDist_ConnectedComponents(LSDIndexRaster& ConnectedComponents);

  ///@details This function takes in a connected components raster and returns an array
  /// of the distance to the baseline for each pixl in the raster
  ///@param ConnectedComponents connected components raster
  ///@return array with baseline components
  ///@author FJC
  ///@date 12/10/17
  Array2D<float> get_DistanceToBaseline_ConnectedComponents(LSDIndexRaster& ConnectedComponents);

  // write profiles to file
  void write_transverse_profile_to_file(LSDRaster& Raster, vector<float> desired_percentiles, float BinWidth, string prefix, int NormaliseToBaseline);
  void write_longitudinal_profile_to_file(LSDRaster& Raster, vector<float> desired_percentiles, float BinWidth, string prefix, int NormaliseToBaseline);
  void print_baseline_to_csv(LSDRaster& ElevationRaster, string csv_filename);
  void write_swath_metadata_to_csv(string csv_filename);
  void write_array_data_to_csv(string csv_filename, LSDFlowInfo& FlowInfo);
  void print_swath_data_to_csvs(string path, string csv_prefix, LSDFlowInfo& FlowInfo, LSDRaster& ElevationRaster);
  // get functions
  // these get data elements
  int get_NPtsInProfile() const {return NPtsInProfile;}
  Array2D<float> get_DistanceToBaselineArray() const { return DistanceToBaselineArray; }
  Array2D<float> get_DistanceAlongBaselineArray() const { return DistanceAlongBaselineArray; }
  Array2D<float> get_BaselineValueArray() const { return BaselineValueArray; }

  float get_XMax() const { return XMax; }
  float get_YMax() const { return YMax; }
  float get_XMin() const { return XMin; }
  float get_YMin() const { return YMin; }
  float get_ProfileHalfWidth() const { return ProfileHalfWidth; }

  vector<int> get_BaselineCols() const { return BaselineCols; }
  vector<int> get_BaselineRows() const { return BaselineRows; }

	protected:

  // Swath template
  vector<float> DistanceAlongBaseline;
  vector<float> BaselineValue;
	vector<int> BaselineRows;  // rows of the baseline points
	vector<int> BaselineCols;  // cols of the baseline points
	Array2D<float> DistanceToBaselineArray;
  Array2D<float> DistanceAlongBaselineArray;
  Array2D<float> BaselineValueArray;

	// metadata
  int NPtsInProfile;
  float ProfileHalfWidth;
  float NoDataValue;
  int NRows;
  int NCols;

  // Bounding Box of profile
  float XMax;
  float XMin;
  float YMax;
  float YMin;

	private:
  void create();
  void create(PointData& ProfilePoints, LSDRaster& RasterTemplate, float ProfileHalfWidth);
  void create(vector<float>& Y_X_points, LSDRaster& RasterTemplate, float& ProfileHalfWidth, float d_space);
  void create(string path, LSDRaster& RasterTemplate, string DEM_prefix, LSDFlowInfo& FlowInfo);



};

#endif
