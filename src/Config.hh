#ifndef __INCLUDE_GUARD_Config_hh__
#define __INCLUDE_GUARD_Config_hh__
#include "Utility.hh"


// Input system:
enum StreamingType { FlatFixed, FlatAdaptive, CurvedFixed, CurvedAdaptive };
inline std::string StreamingName(int n)
{
	std::string name("unknown");
	switch (n)
	{
   		case 0: { name = "FlatFixed";       } break;
   		case 1: { name = "FlatAdaptive";    } break;
   		case 2: { name = "CurvedFixed";     } break;
   		case 3: { name = "CurvedAdaptive";  } break;
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