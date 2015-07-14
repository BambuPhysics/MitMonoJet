#open config file
config = open("config/DMSTreeHeaderConfig.txt", "r")
print config
#print "reading"

# read config file and make lists containing variable names and types
lines = [line for line in config]
variables = []
types = []
tmp = ["name","type"]
for line in lines:
    tmp = line.split()
    if tmp:
        variables.append(tmp[1])
        types.append(tmp[0])
config.close()

# open header file
header = open("Core/interface/MitDMSTree_new.h", "w")
print header
#print "writing"

# write includes and neccesary information
header.write("#ifndef MitDMSTree_H\n")
header.write("#define MitDMSTree_H\n\n")
includes = ["<set>", "<vector>", "<string>", "<utility>", '"TFile.h"', '"TTree.h"', '"TError.h"', '"Math/VectorUtil.h"', '"Math/LorentzVector.h"' ]
for item in includes:
    header.write("#include "+item+"\n"),

# start writing MitDMSTree class
header.write("\nclass MitDMSTree {\n public:\n")
header.write("  /// float doesn't have dictionary by default, so use double\n  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;\n\n")

# Trigger stuff!
header.write("  /// bit map\n  /// DON'T CHANGE ORDER\n  enum Trigger {\n    HLTJetMet   = 1UL<<0,    // event passes jet+met trigger\n    HLTMet      = 1UL<<1,    // event passes met trigger\n    HLTMuon     = 1UL<<2,    // event passes single muon trigger\n    HLTPhoton   = 1UL<<3     // event passes single photon trigger\n  };\n\n")
header.write("  /// bit map\n  /// DON'T CHANGE ORDER\n  enum HLTMatch {\n    JetMatch    = 1UL<<0,    // hardest jet is matched to HLT jet\n    MuonMatch   = 1UL<<1,    // hardest lepton is matched to HLT muon\n    PhotonMatch = 1UL<<2     // hardest photon is matched to HLT photon\n  };\n\n")
header.write("  /// bit map\n  /// DON'T CHANGE ORDER\n  enum Presel {\n    Top      = 1UL<<0,    // event passes top preselection\n    Wlep     = 1UL<<1,    // event passes W>lv preselection\n    Zlep     = 1UL<<2,    // event passes Z>ll preselection\n    Met      = 1UL<<3,    // event passes MET preselection\n    Vbf      = 1UL<<4,    // event passes VBF preselection\n    Gjet     = 1UL<<5     // event passes G+jets preselection\n  };\n\n")

# declare ALL the variables!!!
header.write("  /// variables\n")
for var, type in zip(variables,types):
    declar = "  "+type.ljust(18)+" "+var+"_;\n"
    header.write(declar)

header.write("\n public:\n")
# declare file, tree, and variable vector
header.write("  /// this is the main element\n  TFile *f_;\n  TTree *tree_;\n\n")
header.write("  /// hold the names of variables to facilitate things (filled during Init)\n  std::vector<std::string> variables_;\n\n")

# constructor and deconstructor
header.write("  ///default constructor\n  MitDMSTree():\n    lepPtr1_(&lep1_),lepPtr2_(&lep2_),\n    tauPtr1_(&tau1_),phoPtr1_(&pho1_),\n    fjet1Ptr_(&fjet1_),fjet2Ptr_(&fjet2_),\n    sjetPtr1_(&sjet1_),sjetPtr2_(&sjet2_),\n    jetPtr1_(&jet1_),jetPtr2_(&jet2_),jetPtr3_(&jet3_),jetPtr4_(&jet4_),jetPtr5_(&jet5_),\n    bjetPtr1_(&bjet1_),bjetPtr2_(&bjet2_),\n    genVPtr_(&genV_) {}\n")
header.write("  /// default destructor\n  ~MitDMSTree(){\n    if (f_) f_->Close();\n  };\n\n")

# init vars
header.write("  /// initialize variables and fill list of availables variables\n  void InitVariables();\n\n")

# Load some Trees
header.write('  /// Load a MitDMSTree\n  void LoadTree(const char* file, int type = -1){\n    // to load three different ntutples in the same job DMTree0/1\n    // type == 0/1 if all variables were added\n    // type == -1 (default) if a minimum set of variables was added with tree as name\n    f_ = TFile::Open(file);\n    assert(f_);\n    if     (type == 0) tree_ = dynamic_cast<TTree*>(f_->FindObjectAny("DMSTree"));\n    else               tree_ = dynamic_cast<TTree*>(f_->FindObjectAny("tree"));\n    assert(tree_);\n  }\n\n')
header.write("  /// load a MitDMSTree\n  void LoadTree(const char* file, const char* treeName){\n    f_ = TFile::Open(file);\n    assert(f_);\n    tree_ = dynamic_cast<TTree*>(f_->FindObjectAny(treeName));\n    assert(tree_);\n  }\n\n")

# Start creating a tree
header.write('  /// create a MitDMSTree\n  void CreateTree(int type = -1){\n    assert(type==type); // just to suppress warnings\n    // to create three different ntuples in the same job DMTree0/1\n    // type == 0/1 add all variables\n    // type = -1 (default) add a minimum set of variables with tree as name\n    if     (type == 0) tree_ = new TTree("DMSTree","DM&S ntuples");\n    else               tree_ = new TTree("tree","Smurf ntuples");\n    f_ = 0;\n    InitVariables();\n')

# book the branches of the tree
varType = 'Z'
header.write("    // book the branches\n")
for var,type in zip(variables,types):
    varquote = '"'+var+'",'
    if type == 'LorentzVector':
        varpoint = '&'+var+'Ptr_,'
        vecType = '"ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >", '
        header.write('    tree_->Branch('+varquote.ljust(25)+vecType+varpoint.ljust(25)+');\n')
    else:
        varpoint = '&'+var+'_,'
        if type == "Int_t":
            varType = '"'+var+'/I"'
        elif type == 'UInt_t':
            varType = '"'+var+'/i"'
        elif type == 'float':
            varType = '"'+var+'/F"'
        header.write('    tree_->Branch('+varquote.ljust(25)+varpoint.ljust(25)+varType.ljust(25)+');\n')
header.write("\n  }\n")

# start initializing the tree
header.write("\n  // initialze a MitDMSTree\n  void InitTree(int type = -1){\n    assert(type==type); // just to suppress warnings\n    assert(tree_);\n    // don't forget to set pointers to zero before you set address\n    // or you will fully appreciate that ROOT sucks :)\n    InitVariables();\n    //Set branch address\n    Int_t currentState = gErrorIgnoreLevel;\n    // gErrorIgnoreLevel = kError;\n    gErrorIgnoreLevel = kBreak;\n")

# set some branch addresses
for var,type in zip(variables,types):
    varquote = '"'+var+'",'
    if type == 'LorentzVector':
        varpoint = '&'+var+'Ptr_'
    else:
        varpoint = '&'+var+'_'
    header.write('    tree_->SetBranchAddress('+varquote.ljust(25)+varpoint.ljust(25)+');\n')
header.write('\n    gErrorIgnoreLevel = currentState;\n  }\n\n')

# declare the private lorentz vectors
header.write("\n  private:\n\n")
for var,type in zip(variables,types):
    if type == 'LorentzVector':
        varpoint = '&'+var+'Ptr_'
        header.write("  "+type+"* "+varpoint+";\n")
header.write("\n};\n\n")

# Start initializing the variables
header.write("\ninline void\nMitDMSTree::InitVariables(){\n  // initialize variables\n")
initvalue = '0';
for var, type in zip(variables,types):
    if type == 'UInt_t' or 'Int_t':
        initvalue = '0';
    if type == 'float':
        initvalue = '-999.0'
    if type == 'LorentzVector':
        initvalue = 'LorentzVector()'
    if var == 'puweight':
        initvalue = '1.0'
    var_ = var+"_"
    header.write("  "+var_.ljust(20)+" = "+initvalue.rjust(18)+";\n")

# write end of file!
header.write("}\n\n#endif")
header.close()
