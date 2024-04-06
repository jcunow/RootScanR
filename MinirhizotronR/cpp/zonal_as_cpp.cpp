#include <Rcpp.h>
using namespace Rcpp;

// Custom zonal function
// [[Rcpp::export]]
List custom_zonal_function(RasterLayer raster, RasterLayer zones, Function custom_function) {
  // Get dimensions of the raster and zones
  int nrow = raster.nrow();
  int ncol = raster.ncol();

  // Check if raster dimensions match zone dimensions
  if (nrow != zones.nrow() || ncol != zones.ncol()) {
    stop("Raster and zones dimensions do not match.");
  }

  // Initialize a list to store results for each zone
  std::unordered_map<int, NumericVector> results;

  // Iterate over each cell of the raster
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      double value = raster.get(i, j); // Get raster value at current cell
      int zone = zones.get(i, j); // Get zone value at current cell

      // Skip cells with NA zone values
      if (ISNA(zone))
        continue;

      // Check if zone exists in the results map
      if (results.find(zone) == results.end()) {
        // Create a new NumericVector for the zone if it doesn't exist
        results[zone] = NumericVector::create(value);
      } else {
        // Add value to existing NumericVector for the zone
        results[zone].push_back(value);
      }
    }
  }

  // Apply custom function to each zone vector and store results in a list
  List final_results;
  for (auto it = results.begin(); it != results.end(); ++it) {
    final_results.push_back(custom_function(it->second));
  }

  return final_results;
}
