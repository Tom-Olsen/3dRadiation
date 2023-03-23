# This script generates the code for the template in SphericalHarmonics.cpp:
#
# template<int I>
# double SphericalHarmonicsThPh::Y(double theta, double phi)

def generate_code(L):
    for l in range(L):
        for m in range(-l,l+1):
            i = l*(l+1) + m
            print(f"if constexpr(I == {i:3d}) return Y<{l:>2},{m:>3}>(theta, phi);")
            # print(f"if constexpr(I == {i:3d}) return Y<{l:>2},{m:>3}>(theta, phi);", end="")
        # print("")

generate_code(16)