#include "IO.hpp"

#ifdef CHAINING
std::vector<std::string> col_headers
{
  "cell_type\t",
  "cell_id\t",
  "length\t",
  "radius\t",
  "pos_x\t",  // com position vector components
  "pos_y\t",
  "pos_z\t",
  "ori_x\t",  // orientation vector components
  "ori_y\t",
  "ori_z\t",
  "lower_link\t", // Output who these bacteria are chained to
  "upper_link\n"
};
#else
std::vector<std::string> col_headers
{
  "cell_type\t",
  "cell_id\t",
  "length\t",
  "radius\t",
  "pos_x\t",  // com position vector components
  "pos_y\t",
  "pos_z\t",
  "ori_x\t",  // orientation vector components
  "ori_y\t",
  "ori_z\n"
};
#endif // End chaining

std::string sim_out_dir{
  SIM_DIR
};

void populateCellsFromFile(
  const std::string file_name,
  std::vector< IBacterium* > &element_list,
  bool verbose
)
{
  // Enforce element list to be empty
  element_list.clear();

  std::ifstream inp_file{ file_name, std::ios::in };
  if (!inp_file)
  {
    std::cerr << "Failed to open "
              << file_name
              << " for reading!"
              << std::endl;
    exit (11);
  }
  if ( verbose )
    std::cout << "reading data from " << file_name << '\n';

  // store read lines
  std::string strInput;

  // Print the column titles
  // Should check if the column length is not the same as the correct number
  std::getline(inp_file, strInput);
  std::cout << strInput << '\n';

  std::cout << "Reading in all the same growthFactor" << '\n';
  while (std::getline(inp_file, strInput))
  {
    std::stringstream ss(strInput);

    const long columns { static_cast<long>( col_headers.size() ) };
    std::vector<double> cell_data(col_headers.size());

    std::string cell_type; ss >> cell_type;
    int id; ss>>id;
    double length; ss>>length;
    double radius; ss>>radius;
    Vec3 pos; ss>>pos.x; ss>>pos.y; ss>>pos.z;

    if ( cell_type=="Spherical" )
    {
      std::cout << "Error, not implemented loading spherical yet" << '\n';
      exit(35);
    }
    else if ( cell_type=="RodShaped" )
    {
      std::cout << "attempt to load rod" << '\n';
      Vec3 ori(3); ss>>ori.x; ss>>ori.y; ss>>ori.z;
      double theta { getThetaFromOrientation(ori) };

      element_list.push_back(
        new RodShapedBacterium(
          pos.x,
          pos.y,
          pos.z,
          theta,                              // theta
          0.5*constants::pi,                  // alpha
          constants::nondim_rodGrwthRtePreFac,   // grwthPreFac
          length                        // init_length
        )
      );
    }
    else if ( cell_type=="ChainingRodShaped" )
    {
      std::cout << "attempt to load chaining rod" << '\n';
      Vec3 ori(3); ss>>ori.x; ss>>ori.y; ss>>ori.z;
      double theta { getThetaFromOrientation(ori) };

      element_list.push_back(
        new RodShapedBacterium(
          pos.x,
          pos.y,
          pos.z,
          theta,                              // theta
          0.5*constants::pi,                  // alpha
          constants::nondim_rodGrwthRtePreFac,   // grwthPreFac
          length                        // init_length
        )
      );
      std::cout << "Need to add in springs" << '\n';
      exit(21);
    }
    else
    {
      std::cout << "Cannot load " << cell_type << '\n';
      exit(567);
    }

    exit(23);
  }
  std::cout << "Loaded " << element_list.size() << " cells" << '\n';
}
