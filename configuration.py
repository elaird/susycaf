from supy.defaultConfiguration import *
from supy import whereami

def mainTree() :
    return ("susyTree","tree")

def otherTreesToKeepWhenSkimming() :
    return [("lumiTree","tree")]

def trace() :
    return True

def useCachedFileLists() :
    return True

def cppFiles() :
    return ["cpp/linkdef.cxx"]

def hadd() :
    return ['hadd', whereami()+'/run/phaddy'][1]

def dictionariesToGenerate() :
    return [
        ("pair<string,bool>", "string"),
        ("map<std::string,bool>", "string;map"),
        ("pair<string,string>", "string"),
        ("map<std::string,string>", "string;map"),
        ("ROOT::Math::Cartesian3D<float>", "Math/Point3D.h"),
        ("ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag>", "Math/Vector3D.h"),
        ("vector<ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> >", "vector;Math/Vector3D.h"),
        ("ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag>", "Math/Point3D.h"),
        ("vector<ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> >", "vector;Math/Point3D.h"),
        #ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> > etc. is addressed in linkdef.cxx
        ("vector<ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> > >", "vector;Math/LorentzVector.h"),
        ("vector< vector< float > >", "vector"),
        ("vector< vector< unsigned int> >", "vector"),
        ("vector< vector< int> >", "vector"),
        ]

def detectorSpecs() :
    return {
        "cms": {"etaBE": 1.479, #from CMS PAS EGM-10-005
                "barrelEtaMax": 1.4442,
                "endcapEtaMin": 1.560,
                "CaloSubdetectors": ["Eb", "Ee", "Hbhe", "Hf"],
                "PFSubdetectors": ["Ecal", "Hcal", "Hfem", "Hfhad", "Ps"],
                "CaloRecHitCollections": [""],
                "PFRecHitCollections": ["","cluster"],
                },
        }
