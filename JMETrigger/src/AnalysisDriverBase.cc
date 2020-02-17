#include <Analysis/JMETrigger/interface/AnalysisDriverBase.h>

#include <iostream>

AnalysisDriverBase::AnalysisDriverBase(const std::string& tfile, const std::string& ttree, const std::string& outputFilePath, const std::string& outputFileMode)
  : theFile(tfile.c_str()), theReader(ttree.c_str(), &theFile), outputFilePath_(outputFilePath), outputFileMode_(outputFileMode) {
  map_TTreeReaderValues_.clear();

  auto iter = theReader.GetTree()->GetIteratorOnAllLeaves();
  TLeaf* leaf(nullptr);
  while((leaf = ((TLeaf*) iter->Next()))){
    const std::string type(leaf->GetTypeName());
    if(type == "UInt_t"){
      map_TTreeReaderValues_.insert(std::make_pair(leaf->GetName(), std::make_unique<TTreeReaderValue<unsigned int>>(theReader, leaf->GetName())));
    }
    else if(type == "Int_t"){
      map_TTreeReaderValues_.insert(std::make_pair(leaf->GetName(), std::make_unique<TTreeReaderValue<int>>(theReader, leaf->GetName())));
    }
    else if(type == "Long64_t"){
      map_TTreeReaderValues_.insert(std::make_pair(leaf->GetName(), std::make_unique<TTreeReaderValue<long long>>(theReader, leaf->GetName())));
    }
    else if(type == "ULong64_t"){
      map_TTreeReaderValues_.insert(std::make_pair(leaf->GetName(), std::make_unique<TTreeReaderValue<unsigned long long>>(theReader, leaf->GetName())));
    }
    else if(type == "vector<bool>"){
      map_TTreeReaderValues_.insert(std::make_pair(leaf->GetName(), std::make_unique<TTreeReaderValue<std::vector<bool>>>(theReader, leaf->GetName())));
    }
    else if(type == "vector<unsigned int>"){
      map_TTreeReaderValues_.insert(std::make_pair(leaf->GetName(), std::make_unique<TTreeReaderValue<std::vector<unsigned int>>>(theReader, leaf->GetName())));
    }
    else if(type == "vector<int>"){
      map_TTreeReaderValues_.insert(std::make_pair(leaf->GetName(), std::make_unique<TTreeReaderValue<std::vector<int>>>(theReader, leaf->GetName())));
    }
    else if(type == "vector<float>"){
      map_TTreeReaderValues_.insert(std::make_pair(leaf->GetName(), std::make_unique<TTreeReaderValue<std::vector<float>>>(theReader, leaf->GetName())));
    }
    else if(type == "vector<double>"){
      map_TTreeReaderValues_.insert(std::make_pair(leaf->GetName(), std::make_unique<TTreeReaderValue<std::vector<double>>>(theReader, leaf->GetName())));
    }
    else {
	std::ostringstream ss_str;
      ss_str << "invalid type for input TLeaf \"" << leaf->GetName() << "\": " << type;
      throw std::runtime_error(ss_str.str());
    }
  }
}

AnalysisDriverBase::~AnalysisDriverBase(){
}

bool AnalysisDriverBase::find(const std::string& key) const {

  return (map_TTreeReaderValues_.find(key) != map_TTreeReaderValues_.end());
}

void AnalysisDriverBase::process(const Long64_t firstEntry, const Long64_t maxEntries){

  this->init();

  while(theReader.Next()){

    if(theReader.GetCurrentEntry() < firstEntry){ continue; }

    if(maxEntries == 0){
      break;
    }
    else if(maxEntries > 0){
      if(theReader.GetCurrentEntry() >= (firstEntry+maxEntries)){ continue; }
    }

    this->analyze();
  }

  if(outputFilePath_ != ""){
    this->writeToFile(outputFilePath_, outputFileMode_);
  }
}

void AnalysisDriverBase::writeToFile(const std::string& output_file, const std::string& output_mode){

  TFile outFile(output_file.c_str(), output_mode.c_str());
  if(not outFile.IsZombie()){
    outFile.cd();
    this->write();
  }
  outFile.Close();
}

void AnalysisDriverBase::addOption(const std::string& key, const std::string& opt){

  if(map_options_.find(key) != map_options_.end()){
    std::cout << key << " " << opt << std::endl;
  }
  else {
    map_options_.insert(std::make_pair(key, opt));
  }
}
