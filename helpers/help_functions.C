
int extractClusterId(const std::string& filepath){
	
	size_t pos = filepath.find("cluster");
	if (pos != std::string::npos){

		std::string clusterStr = filepath.substr(pos + 7, 2);
		try {
			return std::stoi(clusterStr);
		} catch (const std::exception& e){
			std::cerr << "Warning: failed to convert cluster number" << clusterStr <<std::endl;
	 	return -1;
		}
	}
	return -1;
}
