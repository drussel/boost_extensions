//  (C) Copyright Daniel Russel 2009. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  copy at http://www.boost.org/LICENSE_1_0.txt)

int main (int argc, char * const argv[]) {
  typedef boost::adjacency_list<boost::vecS, boost::vecS,
                                boost::bidirectionalS,
                                boost::property<boost::vertex_name_t, std::string>,
                                boost::property<boost::edge_weight_t,
                                                double> > Graph;
  Graph pg(10);
  boost::property_map<Graph,
                      boost::vertex_name_t>::type
    vm= boost::get(boost::vertex_name, pg);
  boost::property_map<Graph,
                      boost::edge_weight_t>::type
    ew= boost::get(boost::edge_weight, pg);

  // add edges with weight 1 connecting vertex i with vertex i+5
  for ( int i=0; i< 5; ++i) {
    // sloppy about converting ints to vertex_descriptors
    ew[boost::add_edge(i, i+5, pg).first]=1;
  }
  std::vector<std::pair<int,int> > out
    = get_maximum_weight_bipartite_matching(pg, 5,
                                            boost::get(boost::vertex_index, pg),
                                            boost::get(boost::edge_weight, pg));
  for (unsigned int i=0; i < out.size(); ++i) {
    std::cout << "edge " << out[i].first << "--" << out[i].second << std::endl;
  }
  return 0;
}
