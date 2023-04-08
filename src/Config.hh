#ifndef __INCLUDE_GUARD_Config_hh__
#define __INCLUDE_GUARD_Config_hh__
#include "Utility.hh"


// Input system:
enum StreamingType { FlatStatic, FlatDynamic, CurvedStatic, CurvedDynamic };
inline std::string StreamingName(int n)
{
	std::string name("unknown");
	switch (n)
	{
   		case 0: { name = "FlatStatic";		} break;
   		case 1: { name = "FlatDynamic";		} break;
   		case 2: { name = "CurvedStatic";	} break;
   		case 3: { name = "CurvedDynamic";	} break;
		default: { ExitOnError("Invalid StreamingType"); }
	}
	return name;
}



struct Config
{
    std::string name;
    double simTime;
    int writeFrequency;
    bool updateSphericalHarmonics;
    bool keepSourceNodesActive;
    bool writeData;
    bool printToTerminal;
	bool useCamera;
};
#endif //__INCLUDE_GUARD_Config_hh__