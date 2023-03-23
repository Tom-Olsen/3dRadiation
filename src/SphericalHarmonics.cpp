#include "SphericalHarmonics.h"

// ------------------------------ Spherical Harmonics(theta,phi) ------------------------------
// Associated Legendre Polynomials:
template<int L, int M>
double SphericalHarmonicsThPh::P(double x)
{
    double pmm = 1.0;
    if constexpr(M > 0)
    {
        double somx2 = sqrt((1.0 - x) * (1.0 + x));
        double fact = 1.0;
        for (int i = 1; i <= M; i++)
        {
            pmm *= (-fact) * somx2;
            fact += 2.0;
        }
    }
    if constexpr(L == M)
        return pmm;
    double pmmp1 = x * (2.0 * M + 1.0) * pmm;
    if constexpr(L == M + 1)
        return pmmp1;
    double pll = 0.0;
    for (int ll = M + 2; ll <= L; ++ll)
    {
        pll = ((2.0 * ll - 1.0) * x * pmmp1 - (ll + M - 1.0) * pmm) / (ll - M);
        pmm = pmmp1;
        pmmp1 = pll;
    }
    return pll;
}
// Renormalisation constant for Spherical Harmonic function:
template<int L, int M>
double SphericalHarmonicsThPh::Kd()
{
    double temp = ((2.0 * L + 1.0) * Factorial<L - abs(M)>()) / (4.0 * M_PI * Factorial<L + abs(M)>());
    return sqrt(temp);
}

template<int L, int M>
double SphericalHarmonicsThPh::Y(double theta, double phi)
{
    // l € [0,N]
    // m € [-l,l]
    // theta € [0, pi]
    //   phi € [0,2pi]
    phi = fmod((M * phi) + 200.0 * M_PI, (2.0*M_PI));   // constrain m * phi to [0,2pi]
    constexpr double sqrt2 = 1.414213562373095145;
    if constexpr(M == 0) return Kd<L,0>() * P<L,M>(MyCos(theta));
    else if constexpr(M > 0) return sqrt2 * Kd<L,M>() * MyCos(phi) * P<L,M>(MyCos(theta));
    else return -sqrt2 * Kd<L,-M>() * MySin(phi) * P<L,-M>(MyCos(theta));
}
template<int I>
double SphericalHarmonicsThPh::Y(double theta, double phi)
{
    if constexpr(I ==   0) return Y< 0,  0>(theta, phi);
    if constexpr(I ==   1) return Y< 1, -1>(theta, phi);
    if constexpr(I ==   2) return Y< 1,  0>(theta, phi);
    if constexpr(I ==   3) return Y< 1,  1>(theta, phi);
    if constexpr(I ==   4) return Y< 2, -2>(theta, phi);
    if constexpr(I ==   5) return Y< 2, -1>(theta, phi);
    if constexpr(I ==   6) return Y< 2,  0>(theta, phi);
    if constexpr(I ==   7) return Y< 2,  1>(theta, phi);
    if constexpr(I ==   8) return Y< 2,  2>(theta, phi);
    if constexpr(I ==   9) return Y< 3, -3>(theta, phi);
    if constexpr(I ==  10) return Y< 3, -2>(theta, phi);
    if constexpr(I ==  11) return Y< 3, -1>(theta, phi);
    if constexpr(I ==  12) return Y< 3,  0>(theta, phi);
    if constexpr(I ==  13) return Y< 3,  1>(theta, phi);
    if constexpr(I ==  14) return Y< 3,  2>(theta, phi);
    if constexpr(I ==  15) return Y< 3,  3>(theta, phi);
    if constexpr(I ==  16) return Y< 4, -4>(theta, phi);
    if constexpr(I ==  17) return Y< 4, -3>(theta, phi);
    if constexpr(I ==  18) return Y< 4, -2>(theta, phi);
    if constexpr(I ==  19) return Y< 4, -1>(theta, phi);
    if constexpr(I ==  20) return Y< 4,  0>(theta, phi);
    if constexpr(I ==  21) return Y< 4,  1>(theta, phi);
    if constexpr(I ==  22) return Y< 4,  2>(theta, phi);
    if constexpr(I ==  23) return Y< 4,  3>(theta, phi);
    if constexpr(I ==  24) return Y< 4,  4>(theta, phi);
    if constexpr(I ==  25) return Y< 5, -5>(theta, phi);
    if constexpr(I ==  26) return Y< 5, -4>(theta, phi);
    if constexpr(I ==  27) return Y< 5, -3>(theta, phi);
    if constexpr(I ==  28) return Y< 5, -2>(theta, phi);
    if constexpr(I ==  29) return Y< 5, -1>(theta, phi);
    if constexpr(I ==  30) return Y< 5,  0>(theta, phi);
    if constexpr(I ==  31) return Y< 5,  1>(theta, phi);
    if constexpr(I ==  32) return Y< 5,  2>(theta, phi);
    if constexpr(I ==  33) return Y< 5,  3>(theta, phi);
    if constexpr(I ==  34) return Y< 5,  4>(theta, phi);
    if constexpr(I ==  35) return Y< 5,  5>(theta, phi);
    if constexpr(I ==  36) return Y< 6, -6>(theta, phi);
    if constexpr(I ==  37) return Y< 6, -5>(theta, phi);
    if constexpr(I ==  38) return Y< 6, -4>(theta, phi);
    if constexpr(I ==  39) return Y< 6, -3>(theta, phi);
    if constexpr(I ==  40) return Y< 6, -2>(theta, phi);
    if constexpr(I ==  41) return Y< 6, -1>(theta, phi);
    if constexpr(I ==  42) return Y< 6,  0>(theta, phi);
    if constexpr(I ==  43) return Y< 6,  1>(theta, phi);
    if constexpr(I ==  44) return Y< 6,  2>(theta, phi);
    if constexpr(I ==  45) return Y< 6,  3>(theta, phi);
    if constexpr(I ==  46) return Y< 6,  4>(theta, phi);
    if constexpr(I ==  47) return Y< 6,  5>(theta, phi);
    if constexpr(I ==  48) return Y< 6,  6>(theta, phi);
    if constexpr(I ==  49) return Y< 7, -7>(theta, phi);
    if constexpr(I ==  50) return Y< 7, -6>(theta, phi);
    if constexpr(I ==  51) return Y< 7, -5>(theta, phi);
    if constexpr(I ==  52) return Y< 7, -4>(theta, phi);
    if constexpr(I ==  53) return Y< 7, -3>(theta, phi);
    if constexpr(I ==  54) return Y< 7, -2>(theta, phi);
    if constexpr(I ==  55) return Y< 7, -1>(theta, phi);
    if constexpr(I ==  56) return Y< 7,  0>(theta, phi);
    if constexpr(I ==  57) return Y< 7,  1>(theta, phi);
    if constexpr(I ==  58) return Y< 7,  2>(theta, phi);
    if constexpr(I ==  59) return Y< 7,  3>(theta, phi);
    if constexpr(I ==  60) return Y< 7,  4>(theta, phi);
    if constexpr(I ==  61) return Y< 7,  5>(theta, phi);
    if constexpr(I ==  62) return Y< 7,  6>(theta, phi);
    if constexpr(I ==  63) return Y< 7,  7>(theta, phi);
    if constexpr(I ==  64) return Y< 8, -8>(theta, phi);
    if constexpr(I ==  65) return Y< 8, -7>(theta, phi);
    if constexpr(I ==  66) return Y< 8, -6>(theta, phi);
    if constexpr(I ==  67) return Y< 8, -5>(theta, phi);
    if constexpr(I ==  68) return Y< 8, -4>(theta, phi);
    if constexpr(I ==  69) return Y< 8, -3>(theta, phi);
    if constexpr(I ==  70) return Y< 8, -2>(theta, phi);
    if constexpr(I ==  71) return Y< 8, -1>(theta, phi);
    if constexpr(I ==  72) return Y< 8,  0>(theta, phi);
    if constexpr(I ==  73) return Y< 8,  1>(theta, phi);
    if constexpr(I ==  74) return Y< 8,  2>(theta, phi);
    if constexpr(I ==  75) return Y< 8,  3>(theta, phi);
    if constexpr(I ==  76) return Y< 8,  4>(theta, phi);
    if constexpr(I ==  77) return Y< 8,  5>(theta, phi);
    if constexpr(I ==  78) return Y< 8,  6>(theta, phi);
    if constexpr(I ==  79) return Y< 8,  7>(theta, phi);
    if constexpr(I ==  80) return Y< 8,  8>(theta, phi);
    if constexpr(I ==  81) return Y< 9, -9>(theta, phi);
    if constexpr(I ==  82) return Y< 9, -8>(theta, phi);
    if constexpr(I ==  83) return Y< 9, -7>(theta, phi);
    if constexpr(I ==  84) return Y< 9, -6>(theta, phi);
    if constexpr(I ==  85) return Y< 9, -5>(theta, phi);
    if constexpr(I ==  86) return Y< 9, -4>(theta, phi);
    if constexpr(I ==  87) return Y< 9, -3>(theta, phi);
    if constexpr(I ==  88) return Y< 9, -2>(theta, phi);
    if constexpr(I ==  89) return Y< 9, -1>(theta, phi);
    if constexpr(I ==  90) return Y< 9,  0>(theta, phi);
    if constexpr(I ==  91) return Y< 9,  1>(theta, phi);
    if constexpr(I ==  92) return Y< 9,  2>(theta, phi);
    if constexpr(I ==  93) return Y< 9,  3>(theta, phi);
    if constexpr(I ==  94) return Y< 9,  4>(theta, phi);
    if constexpr(I ==  95) return Y< 9,  5>(theta, phi);
    if constexpr(I ==  96) return Y< 9,  6>(theta, phi);
    if constexpr(I ==  97) return Y< 9,  7>(theta, phi);
    if constexpr(I ==  98) return Y< 9,  8>(theta, phi);
    if constexpr(I ==  99) return Y< 9,  9>(theta, phi);
    if constexpr(I == 100) return Y<10,-10>(theta, phi);
    if constexpr(I == 101) return Y<10, -9>(theta, phi);
    if constexpr(I == 102) return Y<10, -8>(theta, phi);
    if constexpr(I == 103) return Y<10, -7>(theta, phi);
    if constexpr(I == 104) return Y<10, -6>(theta, phi);
    if constexpr(I == 105) return Y<10, -5>(theta, phi);
    if constexpr(I == 106) return Y<10, -4>(theta, phi);
    if constexpr(I == 107) return Y<10, -3>(theta, phi);
    if constexpr(I == 108) return Y<10, -2>(theta, phi);
    if constexpr(I == 109) return Y<10, -1>(theta, phi);
    if constexpr(I == 110) return Y<10,  0>(theta, phi);
    if constexpr(I == 111) return Y<10,  1>(theta, phi);
    if constexpr(I == 112) return Y<10,  2>(theta, phi);
    if constexpr(I == 113) return Y<10,  3>(theta, phi);
    if constexpr(I == 114) return Y<10,  4>(theta, phi);
    if constexpr(I == 115) return Y<10,  5>(theta, phi);
    if constexpr(I == 116) return Y<10,  6>(theta, phi);
    if constexpr(I == 117) return Y<10,  7>(theta, phi);
    if constexpr(I == 118) return Y<10,  8>(theta, phi);
    if constexpr(I == 119) return Y<10,  9>(theta, phi);
    if constexpr(I == 120) return Y<10, 10>(theta, phi);
    if constexpr(I == 121) return Y<11,-11>(theta, phi);
    if constexpr(I == 122) return Y<11,-10>(theta, phi);
    if constexpr(I == 123) return Y<11, -9>(theta, phi);
    if constexpr(I == 124) return Y<11, -8>(theta, phi);
    if constexpr(I == 125) return Y<11, -7>(theta, phi);
    if constexpr(I == 126) return Y<11, -6>(theta, phi);
    if constexpr(I == 127) return Y<11, -5>(theta, phi);
    if constexpr(I == 128) return Y<11, -4>(theta, phi);
    if constexpr(I == 129) return Y<11, -3>(theta, phi);
    if constexpr(I == 130) return Y<11, -2>(theta, phi);
    if constexpr(I == 131) return Y<11, -1>(theta, phi);
    if constexpr(I == 132) return Y<11,  0>(theta, phi);
    if constexpr(I == 133) return Y<11,  1>(theta, phi);
    if constexpr(I == 134) return Y<11,  2>(theta, phi);
    if constexpr(I == 135) return Y<11,  3>(theta, phi);
    if constexpr(I == 136) return Y<11,  4>(theta, phi);
    if constexpr(I == 137) return Y<11,  5>(theta, phi);
    if constexpr(I == 138) return Y<11,  6>(theta, phi);
    if constexpr(I == 139) return Y<11,  7>(theta, phi);
    if constexpr(I == 140) return Y<11,  8>(theta, phi);
    if constexpr(I == 141) return Y<11,  9>(theta, phi);
    if constexpr(I == 142) return Y<11, 10>(theta, phi);
    if constexpr(I == 143) return Y<11, 11>(theta, phi);
    if constexpr(I == 144) return Y<12,-12>(theta, phi);
    if constexpr(I == 145) return Y<12,-11>(theta, phi);
    if constexpr(I == 146) return Y<12,-10>(theta, phi);
    if constexpr(I == 147) return Y<12, -9>(theta, phi);
    if constexpr(I == 148) return Y<12, -8>(theta, phi);
    if constexpr(I == 149) return Y<12, -7>(theta, phi);
    if constexpr(I == 150) return Y<12, -6>(theta, phi);
    if constexpr(I == 151) return Y<12, -5>(theta, phi);
    if constexpr(I == 152) return Y<12, -4>(theta, phi);
    if constexpr(I == 153) return Y<12, -3>(theta, phi);
    if constexpr(I == 154) return Y<12, -2>(theta, phi);
    if constexpr(I == 155) return Y<12, -1>(theta, phi);
    if constexpr(I == 156) return Y<12,  0>(theta, phi);
    if constexpr(I == 157) return Y<12,  1>(theta, phi);
    if constexpr(I == 158) return Y<12,  2>(theta, phi);
    if constexpr(I == 159) return Y<12,  3>(theta, phi);
    if constexpr(I == 160) return Y<12,  4>(theta, phi);
    if constexpr(I == 161) return Y<12,  5>(theta, phi);
    if constexpr(I == 162) return Y<12,  6>(theta, phi);
    if constexpr(I == 163) return Y<12,  7>(theta, phi);
    if constexpr(I == 164) return Y<12,  8>(theta, phi);
    if constexpr(I == 165) return Y<12,  9>(theta, phi);
    if constexpr(I == 166) return Y<12, 10>(theta, phi);
    if constexpr(I == 167) return Y<12, 11>(theta, phi);
    if constexpr(I == 168) return Y<12, 12>(theta, phi);
    if constexpr(I == 169) return Y<13,-13>(theta, phi);
    if constexpr(I == 170) return Y<13,-12>(theta, phi);
    if constexpr(I == 171) return Y<13,-11>(theta, phi);
    if constexpr(I == 172) return Y<13,-10>(theta, phi);
    if constexpr(I == 173) return Y<13, -9>(theta, phi);
    if constexpr(I == 174) return Y<13, -8>(theta, phi);
    if constexpr(I == 175) return Y<13, -7>(theta, phi);
    if constexpr(I == 176) return Y<13, -6>(theta, phi);
    if constexpr(I == 177) return Y<13, -5>(theta, phi);
    if constexpr(I == 178) return Y<13, -4>(theta, phi);
    if constexpr(I == 179) return Y<13, -3>(theta, phi);
    if constexpr(I == 180) return Y<13, -2>(theta, phi);
    if constexpr(I == 181) return Y<13, -1>(theta, phi);
    if constexpr(I == 182) return Y<13,  0>(theta, phi);
    if constexpr(I == 183) return Y<13,  1>(theta, phi);
    if constexpr(I == 184) return Y<13,  2>(theta, phi);
    if constexpr(I == 185) return Y<13,  3>(theta, phi);
    if constexpr(I == 186) return Y<13,  4>(theta, phi);
    if constexpr(I == 187) return Y<13,  5>(theta, phi);
    if constexpr(I == 188) return Y<13,  6>(theta, phi);
    if constexpr(I == 189) return Y<13,  7>(theta, phi);
    if constexpr(I == 190) return Y<13,  8>(theta, phi);
    if constexpr(I == 191) return Y<13,  9>(theta, phi);
    if constexpr(I == 192) return Y<13, 10>(theta, phi);
    if constexpr(I == 193) return Y<13, 11>(theta, phi);
    if constexpr(I == 194) return Y<13, 12>(theta, phi);
    if constexpr(I == 195) return Y<13, 13>(theta, phi);
    if constexpr(I == 196) return Y<14,-14>(theta, phi);
    if constexpr(I == 197) return Y<14,-13>(theta, phi);
    if constexpr(I == 198) return Y<14,-12>(theta, phi);
    if constexpr(I == 199) return Y<14,-11>(theta, phi);
    if constexpr(I == 200) return Y<14,-10>(theta, phi);
    if constexpr(I == 201) return Y<14, -9>(theta, phi);
    if constexpr(I == 202) return Y<14, -8>(theta, phi);
    if constexpr(I == 203) return Y<14, -7>(theta, phi);
    if constexpr(I == 204) return Y<14, -6>(theta, phi);
    if constexpr(I == 205) return Y<14, -5>(theta, phi);
    if constexpr(I == 206) return Y<14, -4>(theta, phi);
    if constexpr(I == 207) return Y<14, -3>(theta, phi);
    if constexpr(I == 208) return Y<14, -2>(theta, phi);
    if constexpr(I == 209) return Y<14, -1>(theta, phi);
    if constexpr(I == 210) return Y<14,  0>(theta, phi);
    if constexpr(I == 211) return Y<14,  1>(theta, phi);
    if constexpr(I == 212) return Y<14,  2>(theta, phi);
    if constexpr(I == 213) return Y<14,  3>(theta, phi);
    if constexpr(I == 214) return Y<14,  4>(theta, phi);
    if constexpr(I == 215) return Y<14,  5>(theta, phi);
    if constexpr(I == 216) return Y<14,  6>(theta, phi);
    if constexpr(I == 217) return Y<14,  7>(theta, phi);
    if constexpr(I == 218) return Y<14,  8>(theta, phi);
    if constexpr(I == 219) return Y<14,  9>(theta, phi);
    if constexpr(I == 220) return Y<14, 10>(theta, phi);
    if constexpr(I == 221) return Y<14, 11>(theta, phi);
    if constexpr(I == 222) return Y<14, 12>(theta, phi);
    if constexpr(I == 223) return Y<14, 13>(theta, phi);
    if constexpr(I == 224) return Y<14, 14>(theta, phi);
    if constexpr(I == 225) return Y<15,-15>(theta, phi);
    if constexpr(I == 226) return Y<15,-14>(theta, phi);
    if constexpr(I == 227) return Y<15,-13>(theta, phi);
    if constexpr(I == 228) return Y<15,-12>(theta, phi);
    if constexpr(I == 229) return Y<15,-11>(theta, phi);
    if constexpr(I == 230) return Y<15,-10>(theta, phi);
    if constexpr(I == 231) return Y<15, -9>(theta, phi);
    if constexpr(I == 232) return Y<15, -8>(theta, phi);
    if constexpr(I == 233) return Y<15, -7>(theta, phi);
    if constexpr(I == 234) return Y<15, -6>(theta, phi);
    if constexpr(I == 235) return Y<15, -5>(theta, phi);
    if constexpr(I == 236) return Y<15, -4>(theta, phi);
    if constexpr(I == 237) return Y<15, -3>(theta, phi);
    if constexpr(I == 238) return Y<15, -2>(theta, phi);
    if constexpr(I == 239) return Y<15, -1>(theta, phi);
    if constexpr(I == 240) return Y<15,  0>(theta, phi);
    if constexpr(I == 241) return Y<15,  1>(theta, phi);
    if constexpr(I == 242) return Y<15,  2>(theta, phi);
    if constexpr(I == 243) return Y<15,  3>(theta, phi);
    if constexpr(I == 244) return Y<15,  4>(theta, phi);
    if constexpr(I == 245) return Y<15,  5>(theta, phi);
    if constexpr(I == 246) return Y<15,  6>(theta, phi);
    if constexpr(I == 247) return Y<15,  7>(theta, phi);
    if constexpr(I == 248) return Y<15,  8>(theta, phi);
    if constexpr(I == 249) return Y<15,  9>(theta, phi);
    if constexpr(I == 250) return Y<15, 10>(theta, phi);
    if constexpr(I == 251) return Y<15, 11>(theta, phi);
    if constexpr(I == 252) return Y<15, 12>(theta, phi);
    if constexpr(I == 253) return Y<15, 13>(theta, phi);
    if constexpr(I == 254) return Y<15, 14>(theta, phi);
    if constexpr(I == 255) return Y<15, 15>(theta, phi);
    if constexpr(I > 255)
    {
        exit_on_error("Spherical Harmonic index i to big.");
        return 0;
    }
}
std::vector<double (*)(double, double)> SphericalHarmonicsThPh::Ylist =
{ &SphericalHarmonicsThPh::Y<  0>, &SphericalHarmonicsThPh::Y<  1>, &SphericalHarmonicsThPh::Y<  2>, &SphericalHarmonicsThPh::Y<  3>, &SphericalHarmonicsThPh::Y<  4>, &SphericalHarmonicsThPh::Y<  5>, &SphericalHarmonicsThPh::Y<  6>, &SphericalHarmonicsThPh::Y<  7>, &SphericalHarmonicsThPh::Y<  8>, &SphericalHarmonicsThPh::Y<  9>,
  &SphericalHarmonicsThPh::Y< 10>, &SphericalHarmonicsThPh::Y< 11>, &SphericalHarmonicsThPh::Y< 12>, &SphericalHarmonicsThPh::Y< 13>, &SphericalHarmonicsThPh::Y< 14>, &SphericalHarmonicsThPh::Y< 15>, &SphericalHarmonicsThPh::Y< 16>, &SphericalHarmonicsThPh::Y< 17>, &SphericalHarmonicsThPh::Y< 18>, &SphericalHarmonicsThPh::Y< 19>,
  &SphericalHarmonicsThPh::Y< 20>, &SphericalHarmonicsThPh::Y< 21>, &SphericalHarmonicsThPh::Y< 22>, &SphericalHarmonicsThPh::Y< 23>, &SphericalHarmonicsThPh::Y< 24>, &SphericalHarmonicsThPh::Y< 25>, &SphericalHarmonicsThPh::Y< 26>, &SphericalHarmonicsThPh::Y< 27>, &SphericalHarmonicsThPh::Y< 28>, &SphericalHarmonicsThPh::Y< 29>,
  &SphericalHarmonicsThPh::Y< 30>, &SphericalHarmonicsThPh::Y< 31>, &SphericalHarmonicsThPh::Y< 32>, &SphericalHarmonicsThPh::Y< 33>, &SphericalHarmonicsThPh::Y< 34>, &SphericalHarmonicsThPh::Y< 35>, &SphericalHarmonicsThPh::Y< 36>, &SphericalHarmonicsThPh::Y< 37>, &SphericalHarmonicsThPh::Y< 38>, &SphericalHarmonicsThPh::Y< 39>,
  &SphericalHarmonicsThPh::Y< 40>, &SphericalHarmonicsThPh::Y< 41>, &SphericalHarmonicsThPh::Y< 42>, &SphericalHarmonicsThPh::Y< 43>, &SphericalHarmonicsThPh::Y< 44>, &SphericalHarmonicsThPh::Y< 45>, &SphericalHarmonicsThPh::Y< 46>, &SphericalHarmonicsThPh::Y< 47>, &SphericalHarmonicsThPh::Y< 48>, &SphericalHarmonicsThPh::Y< 49>,
  &SphericalHarmonicsThPh::Y< 50>, &SphericalHarmonicsThPh::Y< 51>, &SphericalHarmonicsThPh::Y< 52>, &SphericalHarmonicsThPh::Y< 53>, &SphericalHarmonicsThPh::Y< 54>, &SphericalHarmonicsThPh::Y< 55>, &SphericalHarmonicsThPh::Y< 56>, &SphericalHarmonicsThPh::Y< 57>, &SphericalHarmonicsThPh::Y< 58>, &SphericalHarmonicsThPh::Y< 59>,
  &SphericalHarmonicsThPh::Y< 60>, &SphericalHarmonicsThPh::Y< 61>, &SphericalHarmonicsThPh::Y< 62>, &SphericalHarmonicsThPh::Y< 63>, &SphericalHarmonicsThPh::Y< 64>, &SphericalHarmonicsThPh::Y< 65>, &SphericalHarmonicsThPh::Y< 66>, &SphericalHarmonicsThPh::Y< 67>, &SphericalHarmonicsThPh::Y< 68>, &SphericalHarmonicsThPh::Y< 69>,
  &SphericalHarmonicsThPh::Y< 70>, &SphericalHarmonicsThPh::Y< 71>, &SphericalHarmonicsThPh::Y< 72>, &SphericalHarmonicsThPh::Y< 73>, &SphericalHarmonicsThPh::Y< 74>, &SphericalHarmonicsThPh::Y< 75>, &SphericalHarmonicsThPh::Y< 76>, &SphericalHarmonicsThPh::Y< 77>, &SphericalHarmonicsThPh::Y< 78>, &SphericalHarmonicsThPh::Y< 79>,
  &SphericalHarmonicsThPh::Y< 80>, &SphericalHarmonicsThPh::Y< 81>, &SphericalHarmonicsThPh::Y< 82>, &SphericalHarmonicsThPh::Y< 83>, &SphericalHarmonicsThPh::Y< 84>, &SphericalHarmonicsThPh::Y< 85>, &SphericalHarmonicsThPh::Y< 86>, &SphericalHarmonicsThPh::Y< 87>, &SphericalHarmonicsThPh::Y< 88>, &SphericalHarmonicsThPh::Y< 89>,
  &SphericalHarmonicsThPh::Y< 90>, &SphericalHarmonicsThPh::Y< 91>, &SphericalHarmonicsThPh::Y< 92>, &SphericalHarmonicsThPh::Y< 93>, &SphericalHarmonicsThPh::Y< 94>, &SphericalHarmonicsThPh::Y< 95>, &SphericalHarmonicsThPh::Y< 96>, &SphericalHarmonicsThPh::Y< 97>, &SphericalHarmonicsThPh::Y< 98>, &SphericalHarmonicsThPh::Y< 99>,
  &SphericalHarmonicsThPh::Y<100>, &SphericalHarmonicsThPh::Y<101>, &SphericalHarmonicsThPh::Y<102>, &SphericalHarmonicsThPh::Y<103>, &SphericalHarmonicsThPh::Y<104>, &SphericalHarmonicsThPh::Y<105>, &SphericalHarmonicsThPh::Y<106>, &SphericalHarmonicsThPh::Y<107>, &SphericalHarmonicsThPh::Y<108>, &SphericalHarmonicsThPh::Y<109>,
  &SphericalHarmonicsThPh::Y<110>, &SphericalHarmonicsThPh::Y<111>, &SphericalHarmonicsThPh::Y<112>, &SphericalHarmonicsThPh::Y<113>, &SphericalHarmonicsThPh::Y<114>, &SphericalHarmonicsThPh::Y<115>, &SphericalHarmonicsThPh::Y<116>, &SphericalHarmonicsThPh::Y<117>, &SphericalHarmonicsThPh::Y<118>, &SphericalHarmonicsThPh::Y<119>,
  &SphericalHarmonicsThPh::Y<120>, &SphericalHarmonicsThPh::Y<121>, &SphericalHarmonicsThPh::Y<122>, &SphericalHarmonicsThPh::Y<123>, &SphericalHarmonicsThPh::Y<124>, &SphericalHarmonicsThPh::Y<125>, &SphericalHarmonicsThPh::Y<126>, &SphericalHarmonicsThPh::Y<127>, &SphericalHarmonicsThPh::Y<128>, &SphericalHarmonicsThPh::Y<129>,
  &SphericalHarmonicsThPh::Y<130>, &SphericalHarmonicsThPh::Y<131>, &SphericalHarmonicsThPh::Y<132>, &SphericalHarmonicsThPh::Y<133>, &SphericalHarmonicsThPh::Y<134>, &SphericalHarmonicsThPh::Y<135>, &SphericalHarmonicsThPh::Y<136>, &SphericalHarmonicsThPh::Y<137>, &SphericalHarmonicsThPh::Y<138>, &SphericalHarmonicsThPh::Y<139>,
  &SphericalHarmonicsThPh::Y<140>, &SphericalHarmonicsThPh::Y<141>, &SphericalHarmonicsThPh::Y<142>, &SphericalHarmonicsThPh::Y<143>, &SphericalHarmonicsThPh::Y<144>, &SphericalHarmonicsThPh::Y<145>, &SphericalHarmonicsThPh::Y<146>, &SphericalHarmonicsThPh::Y<147>, &SphericalHarmonicsThPh::Y<148>, &SphericalHarmonicsThPh::Y<149>,
  &SphericalHarmonicsThPh::Y<150>, &SphericalHarmonicsThPh::Y<151>, &SphericalHarmonicsThPh::Y<152>, &SphericalHarmonicsThPh::Y<153>, &SphericalHarmonicsThPh::Y<154>, &SphericalHarmonicsThPh::Y<155>, &SphericalHarmonicsThPh::Y<156>, &SphericalHarmonicsThPh::Y<157>, &SphericalHarmonicsThPh::Y<158>, &SphericalHarmonicsThPh::Y<159>,
  &SphericalHarmonicsThPh::Y<160>, &SphericalHarmonicsThPh::Y<161>, &SphericalHarmonicsThPh::Y<162>, &SphericalHarmonicsThPh::Y<163>, &SphericalHarmonicsThPh::Y<164>, &SphericalHarmonicsThPh::Y<165>, &SphericalHarmonicsThPh::Y<166>, &SphericalHarmonicsThPh::Y<167>, &SphericalHarmonicsThPh::Y<168>, &SphericalHarmonicsThPh::Y<169>,
  &SphericalHarmonicsThPh::Y<170>, &SphericalHarmonicsThPh::Y<171>, &SphericalHarmonicsThPh::Y<172>, &SphericalHarmonicsThPh::Y<173>, &SphericalHarmonicsThPh::Y<174>, &SphericalHarmonicsThPh::Y<175>, &SphericalHarmonicsThPh::Y<176>, &SphericalHarmonicsThPh::Y<177>, &SphericalHarmonicsThPh::Y<178>, &SphericalHarmonicsThPh::Y<179>,
  &SphericalHarmonicsThPh::Y<180>, &SphericalHarmonicsThPh::Y<181>, &SphericalHarmonicsThPh::Y<182>, &SphericalHarmonicsThPh::Y<183>, &SphericalHarmonicsThPh::Y<184>, &SphericalHarmonicsThPh::Y<185>, &SphericalHarmonicsThPh::Y<186>, &SphericalHarmonicsThPh::Y<187>, &SphericalHarmonicsThPh::Y<188>, &SphericalHarmonicsThPh::Y<189>,
  &SphericalHarmonicsThPh::Y<190>, &SphericalHarmonicsThPh::Y<191>, &SphericalHarmonicsThPh::Y<192>, &SphericalHarmonicsThPh::Y<193>, &SphericalHarmonicsThPh::Y<194>, &SphericalHarmonicsThPh::Y<195>, &SphericalHarmonicsThPh::Y<196>, &SphericalHarmonicsThPh::Y<197>, &SphericalHarmonicsThPh::Y<198>, &SphericalHarmonicsThPh::Y<199>,
  &SphericalHarmonicsThPh::Y<200>, &SphericalHarmonicsThPh::Y<201>, &SphericalHarmonicsThPh::Y<202>, &SphericalHarmonicsThPh::Y<203>, &SphericalHarmonicsThPh::Y<204>, &SphericalHarmonicsThPh::Y<205>, &SphericalHarmonicsThPh::Y<206>, &SphericalHarmonicsThPh::Y<207>, &SphericalHarmonicsThPh::Y<208>, &SphericalHarmonicsThPh::Y<209>,
  &SphericalHarmonicsThPh::Y<210>, &SphericalHarmonicsThPh::Y<211>, &SphericalHarmonicsThPh::Y<212>, &SphericalHarmonicsThPh::Y<213>, &SphericalHarmonicsThPh::Y<214>, &SphericalHarmonicsThPh::Y<215>, &SphericalHarmonicsThPh::Y<216>, &SphericalHarmonicsThPh::Y<217>, &SphericalHarmonicsThPh::Y<218>, &SphericalHarmonicsThPh::Y<219>,
  &SphericalHarmonicsThPh::Y<220>, &SphericalHarmonicsThPh::Y<221>, &SphericalHarmonicsThPh::Y<222>, &SphericalHarmonicsThPh::Y<223>, &SphericalHarmonicsThPh::Y<224>, &SphericalHarmonicsThPh::Y<225>, &SphericalHarmonicsThPh::Y<226>, &SphericalHarmonicsThPh::Y<227>, &SphericalHarmonicsThPh::Y<228>, &SphericalHarmonicsThPh::Y<229>,
  &SphericalHarmonicsThPh::Y<230>, &SphericalHarmonicsThPh::Y<231>, &SphericalHarmonicsThPh::Y<232>, &SphericalHarmonicsThPh::Y<233>, &SphericalHarmonicsThPh::Y<234>, &SphericalHarmonicsThPh::Y<235>, &SphericalHarmonicsThPh::Y<236>, &SphericalHarmonicsThPh::Y<237>, &SphericalHarmonicsThPh::Y<238>, &SphericalHarmonicsThPh::Y<239>,
  &SphericalHarmonicsThPh::Y<240>, &SphericalHarmonicsThPh::Y<241>, &SphericalHarmonicsThPh::Y<242>, &SphericalHarmonicsThPh::Y<243>, &SphericalHarmonicsThPh::Y<244>, &SphericalHarmonicsThPh::Y<245>, &SphericalHarmonicsThPh::Y<246>, &SphericalHarmonicsThPh::Y<247>, &SphericalHarmonicsThPh::Y<248>, &SphericalHarmonicsThPh::Y<249>,
  &SphericalHarmonicsThPh::Y<250>, &SphericalHarmonicsThPh::Y<251>, &SphericalHarmonicsThPh::Y<252>, &SphericalHarmonicsThPh::Y<253>, &SphericalHarmonicsThPh::Y<254>, &SphericalHarmonicsThPh::Y<255> };



// Y_i(theta,phi), i=l(l+1)+m.
double SphericalHarmonicsThPh::Y(int i, double theta, double phi)
{
    // i = l*(l+1) + m, with l € [0,N], m € [-l,l]
    // l can be reconstructed by taking the sqrt of i and rounding down:    int l = floor(sqrt(i));
    // m can then simply be obtained by rearranged the above equation:      int m = i - l * (l + 1);
    return Ylist[i](theta, phi);
}
// Y_lm(theta,phi).
double SphericalHarmonicsThPh::Y(int l, int m, double theta, double phi)
{
    // l € [0,N]
    // m € [-l,l]
    // theta € [0, pi]
    //   phi € [0,2pi]
    int i = l * (l + 1) + m;
    return Ylist[i](theta, phi);
}



std::vector<double> SphericalHarmonicsThPh::GetCoefficients(const Stencil& stencil, const double* data)
{
    std::vector<double> coefficients(stencil.nCoefficients);
    for(size_t d=0; d<stencil.nDir; d++)
    {
        double theta = stencil.Theta(d);
        double phi = stencil.Phi(d);
        double c = data[d] * stencil.W(d);

        for(size_t i=0; i<stencil.nCoefficients; i++)
            coefficients[i] += c * Y(i,theta,phi);
    }
    return coefficients;
}
void SphericalHarmonicsThPh::GetCoefficients(const Stencil& stencil, const double* data, double* coefficients)
{
    for(size_t i=0; i<stencil.nCoefficients; i++)
        coefficients[i] = 0;
    for(size_t d=0; d<stencil.nDir; d++)
    {
        double theta = stencil.Theta(d);
        double phi = stencil.Phi(d);
        double c = data[d] * stencil.W(d);

        for(size_t i=0; i<stencil.nCoefficients; i++)
            coefficients[i] += c * Y(i,theta,phi);
    }
}
double SphericalHarmonicsThPh::GetValue(double theta, double phi, const std::vector<double>& coefficients, size_t nCoefficients)
{
    double result = 0;
    for(size_t i=0; i<nCoefficients; i++)
        result += coefficients[i] * Y(i,theta,phi);

    return result;
}
double SphericalHarmonicsThPh::GetValue(double theta, double phi, double* coefficients, size_t nCoefficients)
{
    double result = 0;
    for(size_t i=0; i<nCoefficients; i++)
        result += coefficients[i] * Y(i,theta,phi);

    return result;
}
// --------------------------------------------------------------------------------------------



// -------------------------------- Spherical Harmonics(x,y,z) --------------------------------
template<int I>
double SphericalHarmonicsXyz::Y(double x, double y, double z)
{
    // i = l*(l+1) + m, with l € [0,N], m € [-l,l]
    // l can be reconstructed by taking the sqrt of i and rounding down:    int l = floor(sqrt(i));
    // m can then simply be obtained by rearranged the above equation:      int m = i - l * (l + 1);
    constexpr double pi = M_PI;
    if constexpr(I ==   0) return 1/(2.*sqrt(pi));
    if constexpr(I ==   1) return (sqrt(3/pi)*y)/2.;
    if constexpr(I ==   2) return (sqrt(3/pi)*z)/2.;
    if constexpr(I ==   3) return (sqrt(3/pi)*x)/2.;
    if constexpr(I ==   4) return (sqrt(15/pi)*x*y)/2.;
    if constexpr(I ==   5) return (sqrt(15/pi)*y*z)/2.;
    if constexpr(I ==   6) return (sqrt(5/pi)*(-1 + 3*IntegerPow<2>(z)))/4.;
    if constexpr(I ==   7) return (sqrt(15/pi)*x*z)/2.;
    if constexpr(I ==   8) return (sqrt(15/pi)*(x - y)*(x + y))/4.;
    if constexpr(I ==   9) return -(sqrt(35/(2.*pi))*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y)))/4.;
    if constexpr(I ==  10) return (sqrt(105/pi)*x*y*z)/2.;
    if constexpr(I ==  11) return (sqrt(21/(2.*pi))*y*(-1 + 5*IntegerPow<2>(z)))/4.;
    if constexpr(I ==  12) return (sqrt(7/pi)*z*(-3 + 5*IntegerPow<2>(z)))/4.;
    if constexpr(I ==  13) return (sqrt(21/(2.*pi))*x*(-1 + 5*IntegerPow<2>(z)))/4.;
    if constexpr(I ==  14) return (sqrt(105/pi)*(x - y)*(x + y)*z)/4.;
    if constexpr(I ==  15) return (sqrt(35/(2.*pi))*(IntegerPow<3>(x) - 3*x*IntegerPow<2>(y)))/4.;
    if constexpr(I ==  16) return (-3*sqrt(35/pi)*x*y*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z)))/4.;
    if constexpr(I ==  17) return (-3*sqrt(35/(2.*pi))*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*z)/4.;
    if constexpr(I ==  18) return (3*sqrt(5/pi)*x*y*(-1 + 7*IntegerPow<2>(z)))/4.;
    if constexpr(I ==  19) return (3*sqrt(5/(2.*pi))*y*z*(-3 + 7*IntegerPow<2>(z)))/4.;
    if constexpr(I ==  20) return (3*(3 - 30*IntegerPow<2>(z) + 35*IntegerPow<4>(z)))/(16.*sqrt(pi));
    if constexpr(I ==  21) return (3*sqrt(5/(2.*pi))*x*z*(-3 + 7*IntegerPow<2>(z)))/4.;
    if constexpr(I ==  22) return (3*sqrt(5/pi)*(x - y)*(x + y)*(-1 + 7*IntegerPow<2>(z)))/8.;
    if constexpr(I ==  23) return (3*sqrt(35/(2.*pi))*(IntegerPow<3>(x) - 3*x*IntegerPow<2>(y))*z)/4.;
    if constexpr(I ==  24) return (3*sqrt(35/pi)*(IntegerPow<4>(x) - 6*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y)))/16.;
    if constexpr(I ==  25) return (3*sqrt(77/(2.*pi))*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y)))/16.;
    if constexpr(I ==  26) return (-3*sqrt(385/pi)*x*y*z*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z)))/4.;
    if constexpr(I ==  27) return -(sqrt(385/(2.*pi))*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*(-1 + 9*IntegerPow<2>(z)))/16.;
    if constexpr(I ==  28) return (sqrt(1155/pi)*x*y*z*(-1 + 3*IntegerPow<2>(z)))/4.;
    if constexpr(I ==  29) return (sqrt(165/pi)*y*(1 - 14*IntegerPow<2>(z) + 21*IntegerPow<4>(z)))/16.;
    if constexpr(I ==  30) return (sqrt(11/pi)*z*(15 - 70*IntegerPow<2>(z) + 63*IntegerPow<4>(z)))/16.;
    if constexpr(I ==  31) return (sqrt(165/pi)*x*(1 - 14*IntegerPow<2>(z) + 21*IntegerPow<4>(z)))/16.;
    if constexpr(I ==  32) return (sqrt(1155/pi)*(x - y)*(x + y)*z*(-1 + 3*IntegerPow<2>(z)))/8.;
    if constexpr(I ==  33) return (sqrt(385/(2.*pi))*(IntegerPow<3>(x) - 3*x*IntegerPow<2>(y))*(-1 + 9*IntegerPow<2>(z)))/16.;
    if constexpr(I ==  34) return (3*sqrt(385/pi)*(IntegerPow<4>(x) - 6*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*z)/16.;
    if constexpr(I ==  35) return (3*sqrt(77/(2.*pi))*(IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 5*x*IntegerPow<4>(y)))/16.;
    if constexpr(I ==  36) return (sqrt(3003/(2.*pi))*y*(3*IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 3*x*IntegerPow<4>(y)))/16.;
    if constexpr(I ==  37) return (3*sqrt(1001/(2.*pi))*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*z)/16.;
    if constexpr(I ==  38) return (-3*sqrt(91/pi)*x*y*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(-1 + 11*IntegerPow<2>(z)))/8.;
    if constexpr(I ==  39) return -(sqrt(1365/(2.*pi))*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*z*(-3 + 11*IntegerPow<2>(z)))/16.;
    if constexpr(I ==  40) return (sqrt(1365/(2.*pi))*x*y*(1 - 18*IntegerPow<2>(z) + 33*IntegerPow<4>(z)))/16.;
    if constexpr(I ==  41) return (sqrt(273/pi)*y*z*(5 - 30*IntegerPow<2>(z) + 33*IntegerPow<4>(z)))/16.;
    if constexpr(I ==  42) return (sqrt(13/pi)*(-5 + 21*IntegerPow<2>(z)*(5 - 15*IntegerPow<2>(z) + 11*IntegerPow<4>(z))))/32.;
    if constexpr(I ==  43) return (sqrt(273/pi)*x*z*(5 - 30*IntegerPow<2>(z) + 33*IntegerPow<4>(z)))/16.;
    if constexpr(I ==  44) return -(sqrt(1365/(2.*pi))*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(1 - 18*IntegerPow<2>(z) + 33*IntegerPow<4>(z)))/32.;
    if constexpr(I ==  45) return (sqrt(1365/(2.*pi))*(IntegerPow<3>(x) - 3*x*IntegerPow<2>(y))*z*(-3 + 11*IntegerPow<2>(z)))/16.;
    if constexpr(I ==  46) return (3*sqrt(91/pi)*(IntegerPow<4>(x) - 6*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(-1 + 11*IntegerPow<2>(z)))/32.;
    if constexpr(I ==  47) return (3*sqrt(1001/(2.*pi))*(IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 5*x*IntegerPow<4>(y))*z)/16.;
    if constexpr(I ==  48) return (sqrt(3003/(2.*pi))*(IntegerPow<6>(x) - 15*IntegerPow<4>(x)*IntegerPow<2>(y) + 15*IntegerPow<2>(x)*IntegerPow<4>(y) - IntegerPow<6>(y)))/32.;
    if constexpr(I ==  49) return (-3*sqrt(715/pi)*y*(-7*IntegerPow<6>(x) + 35*IntegerPow<4>(x)*IntegerPow<2>(y) - 21*IntegerPow<2>(x)*IntegerPow<4>(y) + IntegerPow<6>(y)))/64.;
    if constexpr(I ==  50) return (3*sqrt(5005/(2.*pi))*y*(3*IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 3*x*IntegerPow<4>(y))*z)/16.;
    if constexpr(I ==  51) return (3*sqrt(385/pi)*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(-1 + 13*IntegerPow<2>(z)))/64.;
    if constexpr(I ==  52) return (-3*sqrt(385/pi)*x*y*z*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(-3 + 13*IntegerPow<2>(z)))/8.;
    if constexpr(I ==  53) return (-3*sqrt(35/pi)*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*(3 - 66*IntegerPow<2>(z) + 143*IntegerPow<4>(z)))/64.;
    if constexpr(I ==  54) return (3*sqrt(35/(2.*pi))*x*y*z*(15 - 110*IntegerPow<2>(z) + 143*IntegerPow<4>(z)))/16.;
    if constexpr(I ==  55) return (sqrt(105/pi)*y*(-5 + 135*IntegerPow<2>(z) - 495*IntegerPow<4>(z) + 429*IntegerPow<6>(z)))/64.;
    if constexpr(I ==  56) return (sqrt(15/pi)*z*(-35 + 315*IntegerPow<2>(z) - 693*IntegerPow<4>(z) + 429*IntegerPow<6>(z)))/32.;
    if constexpr(I ==  57) return (sqrt(105/pi)*x*(-5 + 135*IntegerPow<2>(z) - 495*IntegerPow<4>(z) + 429*IntegerPow<6>(z)))/64.;
    if constexpr(I ==  58) return (-3*sqrt(35/(2.*pi))*z*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(15 - 110*IntegerPow<2>(z) + 143*IntegerPow<4>(z)))/32.;
    if constexpr(I ==  59) return (3*sqrt(35/pi)*(IntegerPow<3>(x) - 3*x*IntegerPow<2>(y))*(3 - 66*IntegerPow<2>(z) + 143*IntegerPow<4>(z)))/64.;
    if constexpr(I ==  60) return (3*sqrt(385/pi)*(IntegerPow<4>(x) - 6*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*z*(-3 + 13*IntegerPow<2>(z)))/32.;
    if constexpr(I ==  61) return (3*sqrt(385/pi)*(IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 5*x*IntegerPow<4>(y))*(-1 + 13*IntegerPow<2>(z)))/64.;
    if constexpr(I ==  62) return (3*sqrt(5005/(2.*pi))*(IntegerPow<6>(x) - 15*IntegerPow<4>(x)*IntegerPow<2>(y) + 15*IntegerPow<2>(x)*IntegerPow<4>(y) - IntegerPow<6>(y))*z)/32.;
    if constexpr(I ==  63) return (3*sqrt(715/pi)*(IntegerPow<7>(x) - 21*IntegerPow<5>(x)*IntegerPow<2>(y) + 35*IntegerPow<3>(x)*IntegerPow<4>(y) - 7*x*IntegerPow<6>(y)))/64.;
    if constexpr(I ==  64) return (-3*sqrt(12155/pi)*y*(-IntegerPow<7>(x) + 7*IntegerPow<5>(x)*IntegerPow<2>(y) - 7*IntegerPow<3>(x)*IntegerPow<4>(y) + x*IntegerPow<6>(y)))/32.;
    if constexpr(I ==  65) return (-3*sqrt(12155/pi)*y*(-7*IntegerPow<6>(x) + 35*IntegerPow<4>(x)*IntegerPow<2>(y) - 21*IntegerPow<2>(x)*IntegerPow<4>(y) + IntegerPow<6>(y))*z)/64.;
    if constexpr(I ==  66) return (sqrt(7293/(2.*pi))*y*(3*IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 3*x*IntegerPow<4>(y))*(-1 + 15*IntegerPow<2>(z)))/32.;
    if constexpr(I ==  67) return (3*sqrt(17017/pi)*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*z*(-1 + 5*IntegerPow<2>(z)))/64.;
    if constexpr(I ==  68) return (-3*sqrt(1309/pi)*x*y*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(1 - 26*IntegerPow<2>(z) + 65*IntegerPow<4>(z)))/32.;
    if constexpr(I ==  69) return -(sqrt(19635/pi)*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*z*(3 - 26*IntegerPow<2>(z) + 39*IntegerPow<4>(z)))/64.;
    if constexpr(I ==  70) return (3*sqrt(595/(2.*pi))*x*y*(-1 + 33*IntegerPow<2>(z) - 143*IntegerPow<4>(z) + 143*IntegerPow<6>(z)))/32.;
    if constexpr(I ==  71) return (3*sqrt(17/pi)*y*z*(-35 + 11*IntegerPow<2>(z)*(35 - 91*IntegerPow<2>(z) + 65*IntegerPow<4>(z))))/64.;
    if constexpr(I ==  72) return (sqrt(17/pi)*(35 - 1260*IntegerPow<2>(z) + 6930*IntegerPow<4>(z) - 12012*IntegerPow<6>(z) + 6435*IntegerPow<8>(z)))/256.;
    if constexpr(I ==  73) return (3*sqrt(17/pi)*x*z*(-35 + 11*IntegerPow<2>(z)*(35 - 91*IntegerPow<2>(z) + 65*IntegerPow<4>(z))))/64.;
    if constexpr(I ==  74) return (3*sqrt(595/(2.*pi))*(x - y)*(x + y)*(-1 + 33*IntegerPow<2>(z) - 143*IntegerPow<4>(z) + 143*IntegerPow<6>(z)))/64.;
    if constexpr(I ==  75) return (sqrt(19635/pi)*(IntegerPow<3>(x) - 3*x*IntegerPow<2>(y))*z*(3 - 26*IntegerPow<2>(z) + 39*IntegerPow<4>(z)))/64.;
    if constexpr(I ==  76) return (3*sqrt(1309/pi)*(1 - 26*IntegerPow<2>(z) + 65*IntegerPow<4>(z))*(8*IntegerPow<4>(y) + 8*IntegerPow<2>(y)*(-1 + IntegerPow<2>(z)) + IntegerPow<2>(-1 + IntegerPow<2>(z))))/128.;
    if constexpr(I ==  77) return (3*sqrt(17017/pi)*(IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 5*x*IntegerPow<4>(y))*z*(-1 + 5*IntegerPow<2>(z)))/64.;
    if constexpr(I ==  78) return (sqrt(7293/(2.*pi))*(IntegerPow<6>(x) - 15*IntegerPow<4>(x)*IntegerPow<2>(y) + 15*IntegerPow<2>(x)*IntegerPow<4>(y) - IntegerPow<6>(y))*(-1 + 15*IntegerPow<2>(z)))/64.;
    if constexpr(I ==  79) return (3*sqrt(12155/pi)*(IntegerPow<7>(x) - 21*IntegerPow<5>(x)*IntegerPow<2>(y) + 35*IntegerPow<3>(x)*IntegerPow<4>(y) - 7*x*IntegerPow<6>(y))*z)/64.;
    if constexpr(I ==  80) return (3*sqrt(12155/pi)*(IntegerPow<8>(x) - 28*IntegerPow<6>(x)*IntegerPow<2>(y) + 70*IntegerPow<4>(x)*IntegerPow<4>(y) - 28*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y)))/256.;
    if constexpr(I ==  81) return (sqrt(230945/(2.*pi))*y*(9*IntegerPow<8>(x) - 84*IntegerPow<6>(x)*IntegerPow<2>(y) + 126*IntegerPow<4>(x)*IntegerPow<4>(y) - 36*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y)))/256.;
    if constexpr(I ==  82) return (-3*sqrt(230945/pi)*y*(-IntegerPow<7>(x) + 7*IntegerPow<5>(x)*IntegerPow<2>(y) - 7*IntegerPow<3>(x)*IntegerPow<4>(y) + x*IntegerPow<6>(y))*z)/32.;
    if constexpr(I ==  83) return (-3*sqrt(13585/(2.*pi))*y*(-7*IntegerPow<6>(x) + 35*IntegerPow<4>(x)*IntegerPow<2>(y) - 21*IntegerPow<2>(x)*IntegerPow<4>(y) + IntegerPow<6>(y))*(-1 + 17*IntegerPow<2>(z)))/256.;
    if constexpr(I ==  84) return (sqrt(40755/(2.*pi))*y*(3*IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 3*x*IntegerPow<4>(y))*z*(-3 + 17*IntegerPow<2>(z)))/32.;
    if constexpr(I ==  85) return (3*sqrt(2717/(2.*pi))*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(1 - 30*IntegerPow<2>(z) + 85*IntegerPow<4>(z)))/128.;
    if constexpr(I ==  86) return (3*sqrt(95095/pi)*x*(x - y)*y*(x + y)*z*(1 - 10*IntegerPow<2>(z) + 17*IntegerPow<4>(z)))/32.;
    if constexpr(I ==  87) return -(sqrt(21945/(2.*pi))*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*(-1 + 39*IntegerPow<2>(z) - 195*IntegerPow<4>(z) + 221*IntegerPow<6>(z)))/128.;
    if constexpr(I ==  88) return (3*sqrt(1045/(2.*pi))*x*y*z*(-7 + 91*IntegerPow<2>(z) - 273*IntegerPow<4>(z) + 221*IntegerPow<6>(z)))/32.;
    if constexpr(I ==  89) return (3*sqrt(95/pi)*y*(7 + 11*IntegerPow<2>(z)*(-28 + 182*IntegerPow<2>(z) - 364*IntegerPow<4>(z) + 221*IntegerPow<6>(z))))/256.;
    if constexpr(I ==  90) return (sqrt(19/pi)*z*(315 - 4620*IntegerPow<2>(z) + 143*IntegerPow<4>(z)*(126 - 180*IntegerPow<2>(z) + 85*IntegerPow<4>(z))))/256.;
    if constexpr(I ==  91) return (3*sqrt(95/pi)*x*(7 + 11*IntegerPow<2>(z)*(-28 + 182*IntegerPow<2>(z) - 364*IntegerPow<4>(z) + 221*IntegerPow<6>(z))))/256.;
    if constexpr(I ==  92) return (3*sqrt(1045/(2.*pi))*(x - y)*(x + y)*z*(-7 + 91*IntegerPow<2>(z) - 273*IntegerPow<4>(z) + 221*IntegerPow<6>(z)))/64.;
    if constexpr(I ==  93) return (sqrt(21945/(2.*pi))*(IntegerPow<3>(x) - 3*x*IntegerPow<2>(y))*(-1 + 39*IntegerPow<2>(z) - 195*IntegerPow<4>(z) + 221*IntegerPow<6>(z)))/128.;
    if constexpr(I ==  94) return (3*sqrt(95095/pi)*z*(1 - 10*IntegerPow<2>(z) + 17*IntegerPow<4>(z))*(8*IntegerPow<4>(y) + 8*IntegerPow<2>(y)*(-1 + IntegerPow<2>(z)) + IntegerPow<2>(-1 + IntegerPow<2>(z))))/128.;
    if constexpr(I ==  95) return (3*sqrt(2717/(2.*pi))*(IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 5*x*IntegerPow<4>(y))*(1 - 30*IntegerPow<2>(z) + 85*IntegerPow<4>(z)))/128.;
    if constexpr(I ==  96) return (sqrt(40755/(2.*pi))*(IntegerPow<6>(x) - 15*IntegerPow<4>(x)*IntegerPow<2>(y) + 15*IntegerPow<2>(x)*IntegerPow<4>(y) - IntegerPow<6>(y))*z*(-3 + 17*IntegerPow<2>(z)))/64.;
    if constexpr(I ==  97) return (3*sqrt(13585/(2.*pi))*(IntegerPow<7>(x) - 21*IntegerPow<5>(x)*IntegerPow<2>(y) + 35*IntegerPow<3>(x)*IntegerPow<4>(y) - 7*x*IntegerPow<6>(y))*(-1 + 17*IntegerPow<2>(z)))/256.;
    if constexpr(I ==  98) return (3*sqrt(230945/pi)*(IntegerPow<8>(x) - 28*IntegerPow<6>(x)*IntegerPow<2>(y) + 70*IntegerPow<4>(x)*IntegerPow<4>(y) - 28*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*z)/256.;
    if constexpr(I ==  99) return (sqrt(230945/(2.*pi))*(IntegerPow<9>(x) - 36*IntegerPow<7>(x)*IntegerPow<2>(y) + 126*IntegerPow<5>(x)*IntegerPow<4>(y) - 84*IntegerPow<3>(x)*IntegerPow<6>(y) + 9*x*IntegerPow<8>(y)))/256.;
    if constexpr(I == 100) return (sqrt(969969/(2.*pi))*x*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + 5*IntegerPow<4>(y)))/256.;
    if constexpr(I == 101) return (sqrt(4849845/(2.*pi))*y*(9*IntegerPow<8>(x) - 84*IntegerPow<6>(x)*IntegerPow<2>(y) + 126*IntegerPow<4>(x)*IntegerPow<4>(y) - 36*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*z)/256.;
    if constexpr(I == 102) return -(sqrt(255255/pi)*y*(-IntegerPow<7>(x) + 7*IntegerPow<5>(x)*IntegerPow<2>(y) - 7*IntegerPow<3>(x)*IntegerPow<4>(y) + x*IntegerPow<6>(y))*(-1 + 19*IntegerPow<2>(z)))/64.;
    if constexpr(I == 103) return (-3*sqrt(85085/(2.*pi))*y*(-7*IntegerPow<6>(x) + 35*IntegerPow<4>(x)*IntegerPow<2>(y) - 21*IntegerPow<2>(x)*IntegerPow<4>(y) + IntegerPow<6>(y))*z*(-3 + 19*IntegerPow<2>(z)))/256.;
    if constexpr(I == 104) return (3*sqrt(5005/(2.*pi))*y*(3*IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 3*x*IntegerPow<4>(y))*(3 - 102*IntegerPow<2>(z) + 323*IntegerPow<4>(z)))/256.;
    if constexpr(I == 105) return (3*sqrt(1001/(2.*pi))*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*z*(15 - 170*IntegerPow<2>(z) + 323*IntegerPow<4>(z)))/128.;
    if constexpr(I == 106) return (-3*sqrt(5005/pi)*x*y*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(-1 + 45*IntegerPow<2>(z) - 255*IntegerPow<4>(z) + 323*IntegerPow<6>(z)))/64.;
    if constexpr(I == 107) return (-3*sqrt(5005/(2.*pi))*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*z*(-7 + 105*IntegerPow<2>(z) - 357*IntegerPow<4>(z) + 323*IntegerPow<6>(z)))/128.;
    if constexpr(I == 108) return (3*sqrt(385/pi)*x*y*(7 + 13*IntegerPow<2>(z)*(-28 + 210*IntegerPow<2>(z) - 476*IntegerPow<4>(z) + 323*IntegerPow<6>(z))))/256.;
    if constexpr(I == 109) return (sqrt(1155/pi)*y*z*(63 + 13*IntegerPow<2>(z)*(-84 + 378*IntegerPow<2>(z) - 612*IntegerPow<4>(z) + 323*IntegerPow<6>(z))))/256.;
    if constexpr(I == 110) return (sqrt(21/pi)*(-63 + 3465*IntegerPow<2>(z) + 143*IntegerPow<4>(z)*(-210 + 630*IntegerPow<2>(z) - 765*IntegerPow<4>(z) + 323*IntegerPow<6>(z))))/512.;
    if constexpr(I == 111) return (sqrt(1155/pi)*x*z*(63 + 13*IntegerPow<2>(z)*(-84 + 378*IntegerPow<2>(z) - 612*IntegerPow<4>(z) + 323*IntegerPow<6>(z))))/256.;
    if constexpr(I == 112) return (-3*sqrt(385/pi)*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(7 + 13*IntegerPow<2>(z)*(-28 + 210*IntegerPow<2>(z) - 476*IntegerPow<4>(z) + 323*IntegerPow<6>(z))))/512.;
    if constexpr(I == 113) return (3*sqrt(5005/(2.*pi))*(IntegerPow<3>(x) - 3*x*IntegerPow<2>(y))*z*(-7 + 105*IntegerPow<2>(z) - 357*IntegerPow<4>(z) + 323*IntegerPow<6>(z)))/128.;
    if constexpr(I == 114) return (3*sqrt(5005/pi)*(IntegerPow<4>(x) - 6*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(-1 + 45*IntegerPow<2>(z) - 255*IntegerPow<4>(z) + 323*IntegerPow<6>(z)))/256.;
    if constexpr(I == 115) return (3*sqrt(1001/(2.*pi))*(IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 5*x*IntegerPow<4>(y))*z*(15 - 170*IntegerPow<2>(z) + 323*IntegerPow<4>(z)))/128.;
    if constexpr(I == 116) return (-3*sqrt(5005/(2.*pi))*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(3 - 102*IntegerPow<2>(z) + 323*IntegerPow<4>(z))*(16*IntegerPow<4>(y) + 16*IntegerPow<2>(y)*(-1 + IntegerPow<2>(z)) + IntegerPow<2>(-1 + IntegerPow<2>(z))))/512.;
    if constexpr(I == 117) return (3*sqrt(85085/(2.*pi))*(IntegerPow<7>(x) - 21*IntegerPow<5>(x)*IntegerPow<2>(y) + 35*IntegerPow<3>(x)*IntegerPow<4>(y) - 7*x*IntegerPow<6>(y))*z*(-3 + 19*IntegerPow<2>(z)))/256.;
    if constexpr(I == 118) return (sqrt(255255/pi)*(IntegerPow<8>(x) - 28*IntegerPow<6>(x)*IntegerPow<2>(y) + 70*IntegerPow<4>(x)*IntegerPow<4>(y) - 28*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*(-1 + 19*IntegerPow<2>(z)))/512.;
    if constexpr(I == 119) return (sqrt(4849845/(2.*pi))*(IntegerPow<9>(x) - 36*IntegerPow<7>(x)*IntegerPow<2>(y) + 126*IntegerPow<5>(x)*IntegerPow<4>(y) - 84*IntegerPow<3>(x)*IntegerPow<6>(y) + 9*x*IntegerPow<8>(y))*z)/256.;
    if constexpr(I == 120) return (sqrt(969969/(2.*pi))*(IntegerPow<10>(x) - 45*IntegerPow<8>(x)*IntegerPow<2>(y) + 210*IntegerPow<6>(x)*IntegerPow<4>(y) - 210*IntegerPow<4>(x)*IntegerPow<6>(y) + 45*IntegerPow<2>(x)*IntegerPow<8>(y) - IntegerPow<10>(y)))/512.;
    if constexpr(I == 121) return -(sqrt(2028117/pi)*y*(-11*IntegerPow<10>(x) + 165*IntegerPow<8>(x)*IntegerPow<2>(y) - 462*IntegerPow<6>(x)*IntegerPow<4>(y) + 330*IntegerPow<4>(x)*IntegerPow<6>(y) - 55*IntegerPow<2>(x)*IntegerPow<8>(y) + IntegerPow<10>(y)))/1024.;
    if constexpr(I == 122) return (sqrt(22309287/(2.*pi))*x*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + 5*IntegerPow<4>(y))*z)/256.;
    if constexpr(I == 123) return (sqrt(1062347/pi)*y*(9*IntegerPow<8>(x) - 84*IntegerPow<6>(x)*IntegerPow<2>(y) + 126*IntegerPow<4>(x)*IntegerPow<4>(y) - 36*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*(-1 + 21*IntegerPow<2>(z)))/1024.;
    if constexpr(I == 124) return -(sqrt(15935205/pi)*y*(-IntegerPow<7>(x) + 7*IntegerPow<5>(x)*IntegerPow<2>(y) - 7*IntegerPow<3>(x)*IntegerPow<4>(y) + x*IntegerPow<6>(y))*z*(-1 + 7*IntegerPow<2>(z)))/64.;
    if constexpr(I == 125) return -(sqrt(838695/pi)*y*(-7*IntegerPow<6>(x) + 35*IntegerPow<4>(x)*IntegerPow<2>(y) - 21*IntegerPow<2>(x)*IntegerPow<4>(y) + IntegerPow<6>(y))*(1 - 38*IntegerPow<2>(z) + 133*IntegerPow<4>(z)))/1024.;
    if constexpr(I == 126) return (sqrt(167739/(2.*pi))*y*(3*IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 3*x*IntegerPow<4>(y))*z*(15 - 190*IntegerPow<2>(z) + 399*IntegerPow<4>(z)))/256.;
    if constexpr(I == 127) return (3*sqrt(3289/pi)*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(-5 + 17*IntegerPow<2>(z)*(15 - 95*IntegerPow<2>(z) + 133*IntegerPow<4>(z))))/1024.;
    if constexpr(I == 128) return (-3*sqrt(23023/pi)*x*y*z*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(-5 + 85*IntegerPow<2>(z) - 323*IntegerPow<4>(z) + 323*IntegerPow<6>(z)))/64.;
    if constexpr(I == 129) return -(sqrt(345345/(2.*pi))*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*(1 - 60*IntegerPow<2>(z) + 510*IntegerPow<4>(z) - 1292*IntegerPow<6>(z) + 969*IntegerPow<8>(z)))/512.;
    if constexpr(I == 130) return (sqrt(49335/pi)*x*y*z*(21 - 420*IntegerPow<2>(z) + 2142*IntegerPow<4>(z) - 3876*IntegerPow<6>(z) + 2261*IntegerPow<8>(z)))/256.;
    if constexpr(I == 131) return (sqrt(759/(2.*pi))*y*(-21 + 13*IntegerPow<2>(z)*(105 - 1050*IntegerPow<2>(z) + 3570*IntegerPow<4>(z) - 4845*IntegerPow<6>(z) + 2261*IntegerPow<8>(z))))/512.;
    if constexpr(I == 132) return (sqrt(23/pi)*z*(-693 + 13*IntegerPow<2>(z)*(1155 - 6930*IntegerPow<2>(z) + 16830*IntegerPow<4>(z) - 17765*IntegerPow<6>(z) + 6783*IntegerPow<8>(z))))/512.;
    if constexpr(I == 133) return (sqrt(759/(2.*pi))*x*(-21 + 13*IntegerPow<2>(z)*(105 - 1050*IntegerPow<2>(z) + 3570*IntegerPow<4>(z) - 4845*IntegerPow<6>(z) + 2261*IntegerPow<8>(z))))/512.;
    if constexpr(I == 134) return -(sqrt(49335/pi)*z*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(21 - 420*IntegerPow<2>(z) + 2142*IntegerPow<4>(z) - 3876*IntegerPow<6>(z) + 2261*IntegerPow<8>(z)))/512.;
    if constexpr(I == 135) return (sqrt(345345/(2.*pi))*(IntegerPow<3>(x) - 3*x*IntegerPow<2>(y))*(1 - 60*IntegerPow<2>(z) + 510*IntegerPow<4>(z) - 1292*IntegerPow<6>(z) + 969*IntegerPow<8>(z)))/512.;
    if constexpr(I == 136) return (3*sqrt(23023/pi)*(IntegerPow<4>(x) - 6*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*z*(-5 + 85*IntegerPow<2>(z) - 323*IntegerPow<4>(z) + 323*IntegerPow<6>(z)))/256.;
    if constexpr(I == 137) return (3*sqrt(3289/pi)*(IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 5*x*IntegerPow<4>(y))*(-5 + 17*IntegerPow<2>(z)*(15 - 95*IntegerPow<2>(z) + 133*IntegerPow<4>(z))))/1024.;
    if constexpr(I == 138) return -(sqrt(167739/(2.*pi))*z*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(15 - 190*IntegerPow<2>(z) + 399*IntegerPow<4>(z))*(16*IntegerPow<4>(y) + 16*IntegerPow<2>(y)*(-1 + IntegerPow<2>(z)) + IntegerPow<2>(-1 + IntegerPow<2>(z))))/512.;
    if constexpr(I == 139) return (sqrt(838695/pi)*(IntegerPow<7>(x) - 21*IntegerPow<5>(x)*IntegerPow<2>(y) + 35*IntegerPow<3>(x)*IntegerPow<4>(y) - 7*x*IntegerPow<6>(y))*(1 - 38*IntegerPow<2>(z) + 133*IntegerPow<4>(z)))/1024.;
    if constexpr(I == 140) return (sqrt(15935205/pi)*(IntegerPow<8>(x) - 28*IntegerPow<6>(x)*IntegerPow<2>(y) + 70*IntegerPow<4>(x)*IntegerPow<4>(y) - 28*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*z*(-1 + 7*IntegerPow<2>(z)))/512.;
    if constexpr(I == 141) return (sqrt(1062347/pi)*(IntegerPow<9>(x) - 36*IntegerPow<7>(x)*IntegerPow<2>(y) + 126*IntegerPow<5>(x)*IntegerPow<4>(y) - 84*IntegerPow<3>(x)*IntegerPow<6>(y) + 9*x*IntegerPow<8>(y))*(-1 + 21*IntegerPow<2>(z)))/1024.;
    if constexpr(I == 142) return (sqrt(22309287/(2.*pi))*(IntegerPow<10>(x) - 45*IntegerPow<8>(x)*IntegerPow<2>(y) + 210*IntegerPow<6>(x)*IntegerPow<4>(y) - 210*IntegerPow<4>(x)*IntegerPow<6>(y) + 45*IntegerPow<2>(x)*IntegerPow<8>(y) - IntegerPow<10>(y))*z)/512.;
    if constexpr(I == 143) return (sqrt(2028117/pi)*(IntegerPow<11>(x) - 55*IntegerPow<9>(x)*IntegerPow<2>(y) + 330*IntegerPow<7>(x)*IntegerPow<4>(y) - 462*IntegerPow<5>(x)*IntegerPow<6>(y) + 165*IntegerPow<3>(x)*IntegerPow<8>(y) - 11*x*IntegerPow<10>(y)))/1024.;
    if constexpr(I == 144) return (5*sqrt(676039/(2.*pi))*x*y*(3*IntegerPow<10>(x) - 55*IntegerPow<8>(x)*IntegerPow<2>(y) + 198*IntegerPow<6>(x)*IntegerPow<4>(y) - 198*IntegerPow<4>(x)*IntegerPow<6>(y) + 55*IntegerPow<2>(x)*IntegerPow<8>(y) - 3*IntegerPow<10>(y)))/512.;
    if constexpr(I == 145) return (-5*sqrt(2028117/pi)*y*(-11*IntegerPow<10>(x) + 165*IntegerPow<8>(x)*IntegerPow<2>(y) - 462*IntegerPow<6>(x)*IntegerPow<4>(y) + 330*IntegerPow<4>(x)*IntegerPow<6>(y) - 55*IntegerPow<2>(x)*IntegerPow<8>(y) + IntegerPow<10>(y))*z)/1024.;
    if constexpr(I == 146) return (5*sqrt(88179/(2.*pi))*x*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + 5*IntegerPow<4>(y))*(-1 + 23*IntegerPow<2>(z)))/512.;
    if constexpr(I == 147) return (5*sqrt(323323/pi)*y*(9*IntegerPow<8>(x) - 84*IntegerPow<6>(x)*IntegerPow<2>(y) + 126*IntegerPow<4>(x)*IntegerPow<4>(y) - 36*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*z*(-3 + 23*IntegerPow<2>(z)))/1024.;
    if constexpr(I == 148) return (-5*sqrt(138567/pi)*y*(-IntegerPow<7>(x) + 7*IntegerPow<5>(x)*IntegerPow<2>(y) - 7*IntegerPow<3>(x)*IntegerPow<4>(y) + x*IntegerPow<6>(y))*(1 - 42*IntegerPow<2>(z) + 161*IntegerPow<4>(z)))/256.;
    if constexpr(I == 149) return (-5*sqrt(138567/pi)*y*(-7*IntegerPow<6>(x) + 35*IntegerPow<4>(x)*IntegerPow<2>(y) - 21*IntegerPow<2>(x)*IntegerPow<4>(y) + IntegerPow<6>(y))*z*(5 - 70*IntegerPow<2>(z) + 161*IntegerPow<4>(z)))/1024.;
    if constexpr(I == 150) return (5*sqrt(2431/(2.*pi))*y*(3*IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 3*x*IntegerPow<4>(y))*(-5 + 285*IntegerPow<2>(z) - 1995*IntegerPow<4>(z) + 3059*IntegerPow<6>(z)))/512.;
    if constexpr(I == 151) return (15*sqrt(17017/pi)*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*z*(-5 + 95*IntegerPow<2>(z) - 399*IntegerPow<4>(z) + 437*IntegerPow<6>(z)))/1024.;
    if constexpr(I == 152) return (-15*sqrt(1001/(2.*pi))*x*y*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(5 + 17*IntegerPow<2>(z)*(-20 + 190*IntegerPow<2>(z) - 532*IntegerPow<4>(z) + 437*IntegerPow<6>(z))))/512.;
    if constexpr(I == 153) return (-5*sqrt(1001/(2.*pi))*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*z*(45 + 17*IntegerPow<2>(z)*(-60 + 342*IntegerPow<2>(z) - 684*IntegerPow<4>(z) + 437*IntegerPow<6>(z))))/512.;
    if constexpr(I == 154) return (5*sqrt(3003/pi)*x*y*(-3 + 225*IntegerPow<2>(z) + 17*IntegerPow<4>(z)*(-150 + 570*IntegerPow<2>(z) - 855*IntegerPow<4>(z) + 437*IntegerPow<6>(z))))/512.;
    if constexpr(I == 155) return (5*sqrt(39/(2.*pi))*y*z*(-231 + 5775*IntegerPow<2>(z) + 17*IntegerPow<4>(z)*(-2310 + 6270*IntegerPow<2>(z) - 7315*IntegerPow<4>(z) + 3059*IntegerPow<6>(z))))/512.;
    if constexpr(I == 156) return (1155 + 65*IntegerPow<2>(z)*(-1386 + 17325*IntegerPow<2>(z) + 17*IntegerPow<4>(z)*(-4620 + 9405*IntegerPow<2>(z) - 8778*IntegerPow<4>(z) + 3059*IntegerPow<6>(z))))/(2048.*sqrt(pi));
    if constexpr(I == 157) return (5*sqrt(39/(2.*pi))*x*z*(-231 + 5775*IntegerPow<2>(z) + 17*IntegerPow<4>(z)*(-2310 + 6270*IntegerPow<2>(z) - 7315*IntegerPow<4>(z) + 3059*IntegerPow<6>(z))))/512.;
    if constexpr(I == 158) return (-5*sqrt(3003/pi)*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(-3 + 225*IntegerPow<2>(z) + 17*IntegerPow<4>(z)*(-150 + 570*IntegerPow<2>(z) - 855*IntegerPow<4>(z) + 437*IntegerPow<6>(z))))/1024.;
    if constexpr(I == 159) return (5*sqrt(1001/(2.*pi))*(IntegerPow<3>(x) - 3*x*IntegerPow<2>(y))*z*(45 + 17*IntegerPow<2>(z)*(-60 + 342*IntegerPow<2>(z) - 684*IntegerPow<4>(z) + 437*IntegerPow<6>(z))))/512.;
    if constexpr(I == 160) return (15*sqrt(1001/(2.*pi))*(8*IntegerPow<4>(y) + 8*IntegerPow<2>(y)*(-1 + IntegerPow<2>(z)) + IntegerPow<2>(-1 + IntegerPow<2>(z)))*(5 + 17*IntegerPow<2>(z)*(-20 + 190*IntegerPow<2>(z) - 532*IntegerPow<4>(z) + 437*IntegerPow<6>(z))))/2048.;
    if constexpr(I == 161) return (15*sqrt(17017/pi)*(IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 5*x*IntegerPow<4>(y))*z*(-5 + 95*IntegerPow<2>(z) - 399*IntegerPow<4>(z) + 437*IntegerPow<6>(z)))/1024.;
    if constexpr(I == 162) return (5*sqrt(2431/(2.*pi))*(IntegerPow<6>(x) - 15*IntegerPow<4>(x)*IntegerPow<2>(y) + 15*IntegerPow<2>(x)*IntegerPow<4>(y) - IntegerPow<6>(y))*(-5 + 285*IntegerPow<2>(z) - 1995*IntegerPow<4>(z) + 3059*IntegerPow<6>(z)))/1024.;
    if constexpr(I == 163) return (5*sqrt(138567/pi)*(IntegerPow<7>(x) - 21*IntegerPow<5>(x)*IntegerPow<2>(y) + 35*IntegerPow<3>(x)*IntegerPow<4>(y) - 7*x*IntegerPow<6>(y))*z*(5 - 70*IntegerPow<2>(z) + 161*IntegerPow<4>(z)))/1024.;
    if constexpr(I == 164) return (5*sqrt(138567/pi)*(1 - 42*IntegerPow<2>(z) + 161*IntegerPow<4>(z))*(128*IntegerPow<8>(y) + 256*IntegerPow<6>(y)*(-1 + IntegerPow<2>(z)) + 160*IntegerPow<4>(y)*IntegerPow<2>(-1 + IntegerPow<2>(z)) + 32*IntegerPow<2>(y)*IntegerPow<3>(-1 + IntegerPow<2>(z)) + IntegerPow<4>(-1 + IntegerPow<2>(z))))/2048.;
    if constexpr(I == 165) return (5*sqrt(323323/pi)*(IntegerPow<9>(x) - 36*IntegerPow<7>(x)*IntegerPow<2>(y) + 126*IntegerPow<5>(x)*IntegerPow<4>(y) - 84*IntegerPow<3>(x)*IntegerPow<6>(y) + 9*x*IntegerPow<8>(y))*z*(-3 + 23*IntegerPow<2>(z)))/1024.;
    if constexpr(I == 166) return (5*sqrt(88179/(2.*pi))*(IntegerPow<10>(x) - 45*IntegerPow<8>(x)*IntegerPow<2>(y) + 210*IntegerPow<6>(x)*IntegerPow<4>(y) - 210*IntegerPow<4>(x)*IntegerPow<6>(y) + 45*IntegerPow<2>(x)*IntegerPow<8>(y) - IntegerPow<10>(y))*(-1 + 23*IntegerPow<2>(z)))/1024.;
    if constexpr(I == 167) return (5*sqrt(2028117/pi)*(IntegerPow<11>(x) - 55*IntegerPow<9>(x)*IntegerPow<2>(y) + 330*IntegerPow<7>(x)*IntegerPow<4>(y) - 462*IntegerPow<5>(x)*IntegerPow<6>(y) + 165*IntegerPow<3>(x)*IntegerPow<8>(y) - 11*x*IntegerPow<10>(y))*z)/1024.;
    if constexpr(I == 168) return (5*sqrt(676039/(2.*pi))*(IntegerPow<12>(x) - 66*IntegerPow<10>(x)*IntegerPow<2>(y) + 495*IntegerPow<8>(x)*IntegerPow<4>(y) - 924*IntegerPow<6>(x)*IntegerPow<6>(y) + 495*IntegerPow<4>(x)*IntegerPow<8>(y) - 66*IntegerPow<2>(x)*IntegerPow<10>(y) + IntegerPow<12>(y)))/2048.;
    if constexpr(I == 169) return (15*sqrt(156009/pi)*y*(13*IntegerPow<12>(x) - 286*IntegerPow<10>(x)*IntegerPow<2>(y) + 1287*IntegerPow<8>(x)*IntegerPow<4>(y) - 1716*IntegerPow<6>(x)*IntegerPow<6>(y) + 715*IntegerPow<4>(x)*IntegerPow<8>(y) - 78*IntegerPow<2>(x)*IntegerPow<10>(y) + IntegerPow<12>(y)))/4096.;
    if constexpr(I == 170) return (15*sqrt(2028117/(2.*pi))*x*y*(3*IntegerPow<10>(x) - 55*IntegerPow<8>(x)*IntegerPow<2>(y) + 198*IntegerPow<6>(x)*IntegerPow<4>(y) - 198*IntegerPow<4>(x)*IntegerPow<6>(y) + 55*IntegerPow<2>(x)*IntegerPow<8>(y) - 3*IntegerPow<10>(y))*z)/512.;
    if constexpr(I == 171) return (-3*sqrt(2028117/pi)*y*(-11*IntegerPow<10>(x) + 165*IntegerPow<8>(x)*IntegerPow<2>(y) - 462*IntegerPow<6>(x)*IntegerPow<4>(y) + 330*IntegerPow<4>(x)*IntegerPow<6>(y) - 55*IntegerPow<2>(x)*IntegerPow<8>(y) + IntegerPow<10>(y))*(-1 + 25*IntegerPow<2>(z)))/4096.;
    if constexpr(I == 172) return (3*sqrt(2028117/(2.*pi))*x*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + 5*IntegerPow<4>(y))*z*(-3 + 25*IntegerPow<2>(z)))/512.;
    if constexpr(I == 173) return (3*sqrt(88179/(2.*pi))*y*(9*IntegerPow<8>(x) - 84*IntegerPow<6>(x)*IntegerPow<2>(y) + 126*IntegerPow<4>(x)*IntegerPow<4>(y) - 36*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*(3 - 138*IntegerPow<2>(z) + 575*IntegerPow<4>(z)))/2048.;
    if constexpr(I == 174) return (-3*sqrt(4849845/pi)*y*(-IntegerPow<7>(x) + 7*IntegerPow<5>(x)*IntegerPow<2>(y) - 7*IntegerPow<3>(x)*IntegerPow<4>(y) + x*IntegerPow<6>(y))*z*(3 - 46*IntegerPow<2>(z) + 115*IntegerPow<4>(z)))/256.;
    if constexpr(I == 175) return (-3*sqrt(692835/(2.*pi))*y*(-7*IntegerPow<6>(x) + 35*IntegerPow<4>(x)*IntegerPow<2>(y) - 21*IntegerPow<2>(x)*IntegerPow<4>(y) + IntegerPow<6>(y))*(-1 + 63*IntegerPow<2>(z) - 483*IntegerPow<4>(z) + 805*IntegerPow<6>(z)))/2048.;
    if constexpr(I == 176) return (3*sqrt(969969/(2.*pi))*y*(3*IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 3*x*IntegerPow<4>(y))*z*(-5 + 105*IntegerPow<2>(z) - 483*IntegerPow<4>(z) + 575*IntegerPow<6>(z)))/512.;
    if constexpr(I == 177) return (3*sqrt(51051/pi)*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(5 + 19*IntegerPow<2>(z)*(-20 + 210*IntegerPow<2>(z) - 644*IntegerPow<4>(z) + 575*IntegerPow<6>(z))))/4096.;
    if constexpr(I == 178) return (-3*sqrt(51051/(2.*pi))*x*y*z*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(45 + 19*IntegerPow<2>(z)*(-60 + 378*IntegerPow<2>(z) - 828*IntegerPow<4>(z) + 575*IntegerPow<6>(z))))/512.;
    if constexpr(I == 179) return (-3*sqrt(15015/pi)*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*(-9 + 765*IntegerPow<2>(z) + 323*IntegerPow<4>(z)*(-30 + 126*IntegerPow<2>(z) - 207*IntegerPow<4>(z) + 115*IntegerPow<6>(z))))/4096.;
    if constexpr(I == 180) return (3*sqrt(1365/pi)*x*y*z*(-99 + 2805*IntegerPow<2>(z) + 323*IntegerPow<4>(z)*(-66 + 198*IntegerPow<2>(z) - 253*IntegerPow<4>(z) + 115*IntegerPow<6>(z))))/512.;
    if constexpr(I == 181) return (3*sqrt(273/pi)*y*(33 - 2970*IntegerPow<2>(z) + 42075*IntegerPow<4>(z) + 323*IntegerPow<6>(z)*(-660 + 1485*IntegerPow<2>(z) - 1518*IntegerPow<4>(z) + 575*IntegerPow<6>(z))))/2048.;
    if constexpr(I == 182) return (3*sqrt(3/pi)*z*(3003 - 90090*IntegerPow<2>(z) + 765765*IntegerPow<4>(z) + 323*IntegerPow<6>(z)*(-8580 + 7*IntegerPow<2>(z)*(2145 - 1794*IntegerPow<2>(z) + 575*IntegerPow<4>(z)))))/2048.;
    if constexpr(I == 183) return (3*sqrt(273/pi)*x*(33 - 2970*IntegerPow<2>(z) + 42075*IntegerPow<4>(z) + 323*IntegerPow<6>(z)*(-660 + 1485*IntegerPow<2>(z) - 1518*IntegerPow<4>(z) + 575*IntegerPow<6>(z))))/2048.;
    if constexpr(I == 184) return (-3*sqrt(1365/pi)*z*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(-99 + 2805*IntegerPow<2>(z) + 323*IntegerPow<4>(z)*(-66 + 198*IntegerPow<2>(z) - 253*IntegerPow<4>(z) + 115*IntegerPow<6>(z))))/1024.;
    if constexpr(I == 185) return (3*sqrt(15015/pi)*(IntegerPow<3>(x) - 3*x*IntegerPow<2>(y))*(-9 + 765*IntegerPow<2>(z) + 323*IntegerPow<4>(z)*(-30 + 126*IntegerPow<2>(z) - 207*IntegerPow<4>(z) + 115*IntegerPow<6>(z))))/4096.;
    if constexpr(I == 186) return (3*sqrt(51051/(2.*pi))*z*(8*IntegerPow<4>(y) + 8*IntegerPow<2>(y)*(-1 + IntegerPow<2>(z)) + IntegerPow<2>(-1 + IntegerPow<2>(z)))*(45 + 19*IntegerPow<2>(z)*(-60 + 378*IntegerPow<2>(z) - 828*IntegerPow<4>(z) + 575*IntegerPow<6>(z))))/2048.;
    if constexpr(I == 187) return (3*sqrt(51051/pi)*(IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 5*x*IntegerPow<4>(y))*(5 + 19*IntegerPow<2>(z)*(-20 + 210*IntegerPow<2>(z) - 644*IntegerPow<4>(z) + 575*IntegerPow<6>(z))))/4096.;
    if constexpr(I == 188) return (3*sqrt(969969/(2.*pi))*(IntegerPow<6>(x) - 15*IntegerPow<4>(x)*IntegerPow<2>(y) + 15*IntegerPow<2>(x)*IntegerPow<4>(y) - IntegerPow<6>(y))*z*(-5 + 105*IntegerPow<2>(z) - 483*IntegerPow<4>(z) + 575*IntegerPow<6>(z)))/1024.;
    if constexpr(I == 189) return (3*sqrt(692835/(2.*pi))*(IntegerPow<7>(x) - 21*IntegerPow<5>(x)*IntegerPow<2>(y) + 35*IntegerPow<3>(x)*IntegerPow<4>(y) - 7*x*IntegerPow<6>(y))*(-1 + 63*IntegerPow<2>(z) - 483*IntegerPow<4>(z) + 805*IntegerPow<6>(z)))/2048.;
    if constexpr(I == 190) return (3*sqrt(4849845/pi)*z*(3 - 46*IntegerPow<2>(z) + 115*IntegerPow<4>(z))*(128*IntegerPow<8>(y) + 256*IntegerPow<6>(y)*(-1 + IntegerPow<2>(z)) + 160*IntegerPow<4>(y)*IntegerPow<2>(-1 + IntegerPow<2>(z)) + 32*IntegerPow<2>(y)*IntegerPow<3>(-1 + IntegerPow<2>(z)) + IntegerPow<4>(-1 + IntegerPow<2>(z))))/2048.;
    if constexpr(I == 191) return (3*sqrt(88179/(2.*pi))*(IntegerPow<9>(x) - 36*IntegerPow<7>(x)*IntegerPow<2>(y) + 126*IntegerPow<5>(x)*IntegerPow<4>(y) - 84*IntegerPow<3>(x)*IntegerPow<6>(y) + 9*x*IntegerPow<8>(y))*(3 - 138*IntegerPow<2>(z) + 575*IntegerPow<4>(z)))/2048.;
    if constexpr(I == 192) return (3*sqrt(2028117/(2.*pi))*(IntegerPow<10>(x) - 45*IntegerPow<8>(x)*IntegerPow<2>(y) + 210*IntegerPow<6>(x)*IntegerPow<4>(y) - 210*IntegerPow<4>(x)*IntegerPow<6>(y) + 45*IntegerPow<2>(x)*IntegerPow<8>(y) - IntegerPow<10>(y))*z*(-3 + 25*IntegerPow<2>(z)))/1024.;
    if constexpr(I == 193) return (3*sqrt(2028117/pi)*(IntegerPow<11>(x) - 55*IntegerPow<9>(x)*IntegerPow<2>(y) + 330*IntegerPow<7>(x)*IntegerPow<4>(y) - 462*IntegerPow<5>(x)*IntegerPow<6>(y) + 165*IntegerPow<3>(x)*IntegerPow<8>(y) - 11*x*IntegerPow<10>(y))*(-1 + 25*IntegerPow<2>(z)))/4096.;
    if constexpr(I == 194) return (15*sqrt(2028117/(2.*pi))*(IntegerPow<12>(x) - 66*IntegerPow<10>(x)*IntegerPow<2>(y) + 495*IntegerPow<8>(x)*IntegerPow<4>(y) - 924*IntegerPow<6>(x)*IntegerPow<6>(y) + 495*IntegerPow<4>(x)*IntegerPow<8>(y) - 66*IntegerPow<2>(x)*IntegerPow<10>(y) + IntegerPow<12>(y))*z)/2048.;
    if constexpr(I == 195) return (15*sqrt(156009/pi)*(IntegerPow<13>(x) - 78*IntegerPow<11>(x)*IntegerPow<2>(y) + 715*IntegerPow<9>(x)*IntegerPow<4>(y) - 1716*IntegerPow<7>(x)*IntegerPow<6>(y) + 1287*IntegerPow<5>(x)*IntegerPow<8>(y) - 286*IntegerPow<3>(x)*IntegerPow<10>(y) + 13*x*IntegerPow<12>(y)))/4096.;
    if constexpr(I == 196) return (15*sqrt(646323/pi)*x*y*(IntegerPow<6>(x) - 21*IntegerPow<4>(x)*IntegerPow<2>(y) + 35*IntegerPow<2>(x)*IntegerPow<4>(y) - 7*IntegerPow<6>(y))*(7*IntegerPow<6>(x) - 35*IntegerPow<4>(x)*IntegerPow<2>(y) + 21*IntegerPow<2>(x)*IntegerPow<4>(y) - IntegerPow<6>(y)))/4096.;
    if constexpr(I == 197) return (15*sqrt(4524261/pi)*y*(13*IntegerPow<12>(x) - 286*IntegerPow<10>(x)*IntegerPow<2>(y) + 1287*IntegerPow<8>(x)*IntegerPow<4>(y) - 1716*IntegerPow<6>(x)*IntegerPow<6>(y) + 715*IntegerPow<4>(x)*IntegerPow<8>(y) - 78*IntegerPow<2>(x)*IntegerPow<10>(y) + IntegerPow<12>(y))*z)/4096.;
    if constexpr(I == 198) return (5*sqrt(1508087/(2.*pi))*x*(x - y)*y*(x + y)*(IntegerPow<2>(x) - 3*IntegerPow<2>(y))*(3*IntegerPow<2>(x) - IntegerPow<2>(y))*(IntegerPow<2>(x) - 4*x*y + IntegerPow<2>(y))*(IntegerPow<2>(x) + 4*x*y + IntegerPow<2>(y))*(-1 + 27*IntegerPow<2>(z)))/1024.;
    if constexpr(I == 199) return (-5*sqrt(58815393/pi)*y*(-11*IntegerPow<10>(x) + 165*IntegerPow<8>(x)*IntegerPow<2>(y) - 462*IntegerPow<6>(x)*IntegerPow<4>(y) + 330*IntegerPow<4>(x)*IntegerPow<6>(y) - 55*IntegerPow<2>(x)*IntegerPow<8>(y) + IntegerPow<10>(y))*z*(-1 + 9*IntegerPow<2>(z)))/4096.;
    if constexpr(I == 200) return (sqrt(58815393/pi)*x*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + 5*IntegerPow<4>(y))*(-1 + 5*IntegerPow<2>(z))*(-1 + 45*IntegerPow<2>(z)))/4096.;
    if constexpr(I == 201) return (sqrt(98025655/(2.*pi))*y*(9*IntegerPow<8>(x) - 84*IntegerPow<6>(x)*IntegerPow<2>(y) + 126*IntegerPow<4>(x)*IntegerPow<4>(y) - 36*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*z*(3 - 50*IntegerPow<2>(z) + 135*IntegerPow<4>(z)))/2048.;
    if constexpr(I == 202) return -(sqrt(12785955/pi)*y*(-IntegerPow<7>(x) + 7*IntegerPow<5>(x)*IntegerPow<2>(y) - 7*IntegerPow<3>(x)*IntegerPow<4>(y) + x*IntegerPow<6>(y))*(-1 + 23*IntegerPow<2>(z)*(3 - 25*IntegerPow<2>(z) + 45*IntegerPow<4>(z))))/512.;
    if constexpr(I == 203) return -(sqrt(20092215/(2.*pi))*y*(-7*IntegerPow<6>(x) + 35*IntegerPow<4>(x)*IntegerPow<2>(y) - 21*IntegerPow<2>(x)*IntegerPow<4>(y) + IntegerPow<6>(y))*z*(-7 + 23*IntegerPow<2>(z)*(7 - 35*IntegerPow<2>(z) + 45*IntegerPow<4>(z))))/2048.;
    if constexpr(I == 204) return (sqrt(46881835/pi)*y*(3*IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 3*x*IntegerPow<4>(y))*(1 - 84*IntegerPow<2>(z) + 966*IntegerPow<4>(z) - 3220*IntegerPow<6>(z) + 3105*IntegerPow<8>(z)))/4096.;
    if constexpr(I == 205) return (3*sqrt(9376367/pi)*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*z*(-1 + 5*IntegerPow<2>(z))*(-5 + 23*IntegerPow<2>(z)*(5 - 17*IntegerPow<2>(z) + 15*IntegerPow<4>(z))))/4096.;
    if constexpr(I == 206) return (-3*sqrt(2467465/(2.*pi))*x*y*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(-1 + 19*IntegerPow<2>(z)*(5 - 70*IntegerPow<2>(z) + 322*IntegerPow<4>(z) - 575*IntegerPow<6>(z) + 345*IntegerPow<8>(z))))/1024.;
    if constexpr(I == 207) return -(sqrt(224315/pi)*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*z*(-99 + 19*IntegerPow<2>(z)*(165 - 1386*IntegerPow<2>(z) + 4554*IntegerPow<4>(z) - 6325*IntegerPow<6>(z) + 3105*IntegerPow<8>(z))))/4096.;
    if constexpr(I == 208) return (sqrt(39585/pi)*x*y*(33 - 3366*IntegerPow<2>(z) + 323*IntegerPow<4>(z)*(165 - 924*IntegerPow<2>(z) + 23*IntegerPow<4>(z)*(99 - 110*IntegerPow<2>(z) + 45*IntegerPow<4>(z)))))/4096.;
    if constexpr(I == 209) return (sqrt(3045/pi)*y*z*(429 + 17*IntegerPow<2>(z)*(-858 + 19*IntegerPow<2>(z)*(429 - 1716*IntegerPow<2>(z) + 3289*IntegerPow<4>(z) - 2990*IntegerPow<6>(z) + 1035*IntegerPow<8>(z)))))/2048.;
    if constexpr(I == 210) return (sqrt(29/pi)*(-429 + 45045*IntegerPow<2>(z) - 765765*IntegerPow<4>(z) + 323*IntegerPow<6>(z)*(15015 - 45045*IntegerPow<2>(z) + 69069*IntegerPow<4>(z) - 52325*IntegerPow<6>(z) + 15525*IntegerPow<8>(z))))/4096.;
    if constexpr(I == 211) return (sqrt(3045/pi)*x*z*(429 + 17*IntegerPow<2>(z)*(-858 + 19*IntegerPow<2>(z)*(429 - 1716*IntegerPow<2>(z) + 3289*IntegerPow<4>(z) - 2990*IntegerPow<6>(z) + 1035*IntegerPow<8>(z)))))/2048.;
    if constexpr(I == 212) return -(sqrt(39585/pi)*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(33 - 3366*IntegerPow<2>(z) + 323*IntegerPow<4>(z)*(165 - 924*IntegerPow<2>(z) + 23*IntegerPow<4>(z)*(99 - 110*IntegerPow<2>(z) + 45*IntegerPow<4>(z)))))/8192.;
    if constexpr(I == 213) return (sqrt(224315/pi)*(IntegerPow<3>(x) - 3*x*IntegerPow<2>(y))*z*(-99 + 19*IntegerPow<2>(z)*(165 - 1386*IntegerPow<2>(z) + 4554*IntegerPow<4>(z) - 6325*IntegerPow<6>(z) + 3105*IntegerPow<8>(z))))/4096.;
    if constexpr(I == 214) return (3*sqrt(2467465/(2.*pi))*(8*IntegerPow<4>(y) + 8*IntegerPow<2>(y)*(-1 + IntegerPow<2>(z)) + IntegerPow<2>(-1 + IntegerPow<2>(z)))*(-1 + 19*IntegerPow<2>(z)*(5 - 70*IntegerPow<2>(z) + 322*IntegerPow<4>(z) - 575*IntegerPow<6>(z) + 345*IntegerPow<8>(z))))/4096.;
    if constexpr(I == 215) return (3*sqrt(9376367/pi)*x*(IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + 5*IntegerPow<4>(y))*z*(-1 + 5*IntegerPow<2>(z))*(-5 + 23*IntegerPow<2>(z)*(5 - 17*IntegerPow<2>(z) + 15*IntegerPow<4>(z))))/4096.;
    if constexpr(I == 216) return -(sqrt(46881835/pi)*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(1 - 84*IntegerPow<2>(z) + 966*IntegerPow<4>(z) - 3220*IntegerPow<6>(z) + 3105*IntegerPow<8>(z))*(16*IntegerPow<4>(y) + 16*IntegerPow<2>(y)*(-1 + IntegerPow<2>(z)) + IntegerPow<2>(-1 + IntegerPow<2>(z))))/8192.;
    if constexpr(I == 217) return (sqrt(20092215/(2.*pi))*(IntegerPow<7>(x) - 21*IntegerPow<5>(x)*IntegerPow<2>(y) + 35*IntegerPow<3>(x)*IntegerPow<4>(y) - 7*x*IntegerPow<6>(y))*z*(-7 + 23*IntegerPow<2>(z)*(7 - 35*IntegerPow<2>(z) + 45*IntegerPow<4>(z))))/2048.;
    if constexpr(I == 218) return (sqrt(12785955/pi)*(IntegerPow<8>(x) - 28*IntegerPow<6>(x)*IntegerPow<2>(y) + 70*IntegerPow<4>(x)*IntegerPow<4>(y) - 28*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*(-1 + 23*IntegerPow<2>(z)*(3 - 25*IntegerPow<2>(z) + 45*IntegerPow<4>(z))))/4096.;
    if constexpr(I == 219) return (sqrt(98025655/(2.*pi))*(IntegerPow<9>(x) - 36*IntegerPow<7>(x)*IntegerPow<2>(y) + 126*IntegerPow<5>(x)*IntegerPow<4>(y) - 84*IntegerPow<3>(x)*IntegerPow<6>(y) + 9*x*IntegerPow<8>(y))*z*(3 - 50*IntegerPow<2>(z) + 135*IntegerPow<4>(z)))/2048.;
    if constexpr(I == 220) return -(sqrt(58815393/pi)*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(-1 + 5*IntegerPow<2>(z))*(-1 + 45*IntegerPow<2>(z))*(256*IntegerPow<8>(y) + 512*IntegerPow<6>(y)*(-1 + IntegerPow<2>(z)) + 304*IntegerPow<4>(y)*IntegerPow<2>(-1 + IntegerPow<2>(z)) + 48*IntegerPow<2>(y)*IntegerPow<3>(-1 + IntegerPow<2>(z)) + IntegerPow<4>(-1 + IntegerPow<2>(z))))/8192.;
    if constexpr(I == 221) return (5*sqrt(58815393/pi)*(IntegerPow<11>(x) - 55*IntegerPow<9>(x)*IntegerPow<2>(y) + 330*IntegerPow<7>(x)*IntegerPow<4>(y) - 462*IntegerPow<5>(x)*IntegerPow<6>(y) + 165*IntegerPow<3>(x)*IntegerPow<8>(y) - 11*x*IntegerPow<10>(y))*z*(-1 + 9*IntegerPow<2>(z)))/4096.;
    if constexpr(I == 222) return (5*sqrt(1508087/(2.*pi))*(IntegerPow<12>(x) - 66*IntegerPow<10>(x)*IntegerPow<2>(y) + 495*IntegerPow<8>(x)*IntegerPow<4>(y) - 924*IntegerPow<6>(x)*IntegerPow<6>(y) + 495*IntegerPow<4>(x)*IntegerPow<8>(y) - 66*IntegerPow<2>(x)*IntegerPow<10>(y) + IntegerPow<12>(y))*(-1 + 27*IntegerPow<2>(z)))/4096.;
    if constexpr(I == 223) return (15*sqrt(4524261/pi)*(IntegerPow<13>(x) - 78*IntegerPow<11>(x)*IntegerPow<2>(y) + 715*IntegerPow<9>(x)*IntegerPow<4>(y) - 1716*IntegerPow<7>(x)*IntegerPow<6>(y) + 1287*IntegerPow<5>(x)*IntegerPow<8>(y) - 286*IntegerPow<3>(x)*IntegerPow<10>(y) + 13*x*IntegerPow<12>(y))*z)/4096.;
    if constexpr(I == 224) return (15*sqrt(646323/pi)*(IntegerPow<14>(x) - 91*IntegerPow<12>(x)*IntegerPow<2>(y) + 1001*IntegerPow<10>(x)*IntegerPow<4>(y) - 3003*IntegerPow<8>(x)*IntegerPow<6>(y) + 3003*IntegerPow<6>(x)*IntegerPow<8>(y) - 1001*IntegerPow<4>(x)*IntegerPow<10>(y) + 91*IntegerPow<2>(x)*IntegerPow<12>(y) - IntegerPow<14>(y)))/8192.;
    if constexpr(I == 225) return (-3*sqrt(33393355/(2.*pi))*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(IntegerPow<8>(x) - 28*IntegerPow<6>(x)*IntegerPow<2>(y) + 134*IntegerPow<4>(x)*IntegerPow<4>(y) - 92*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y)))/8192.;
    if constexpr(I == 226) return (15*sqrt(20036013/pi)*x*y*(IntegerPow<6>(x) - 21*IntegerPow<4>(x)*IntegerPow<2>(y) + 35*IntegerPow<2>(x)*IntegerPow<4>(y) - 7*IntegerPow<6>(y))*(7*IntegerPow<6>(x) - 35*IntegerPow<4>(x)*IntegerPow<2>(y) + 21*IntegerPow<2>(x)*IntegerPow<4>(y) - IntegerPow<6>(y))*z)/4096.;
    if constexpr(I == 227) return (15*sqrt(690897/(2.*pi))*y*(13*IntegerPow<12>(x) - 286*IntegerPow<10>(x)*IntegerPow<2>(y) + 1287*IntegerPow<8>(x)*IntegerPow<4>(y) - 1716*IntegerPow<6>(x)*IntegerPow<6>(y) + 715*IntegerPow<4>(x)*IntegerPow<8>(y) - 78*IntegerPow<2>(x)*IntegerPow<10>(y) + IntegerPow<12>(y))*(-1 + 29*IntegerPow<2>(z)))/8192.;
    if constexpr(I == 228) return (15*sqrt(1612093/(2.*pi))*x*(x - y)*y*(x + y)*(IntegerPow<2>(x) - 3*IntegerPow<2>(y))*(3*IntegerPow<2>(x) - IntegerPow<2>(y))*(IntegerPow<2>(x) - 4*x*y + IntegerPow<2>(y))*(IntegerPow<2>(x) + 4*x*y + IntegerPow<2>(y))*z*(-3 + 29*IntegerPow<2>(z)))/1024.;
    if constexpr(I == 229) return (-5*sqrt(4836279/(2.*pi))*y*(-11*IntegerPow<10>(x) + 165*IntegerPow<8>(x)*IntegerPow<2>(y) - 462*IntegerPow<6>(x)*IntegerPow<4>(y) + 330*IntegerPow<4>(x)*IntegerPow<6>(y) - 55*IntegerPow<2>(x)*IntegerPow<8>(y) + IntegerPow<10>(y))*(1 - 54*IntegerPow<2>(z) + 261*IntegerPow<4>(z)))/8192.;
    if constexpr(I == 230) return (sqrt(314358135/pi)*x*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + 5*IntegerPow<4>(y))*z*(5 - 90*IntegerPow<2>(z) + 261*IntegerPow<4>(z)))/4096.;
    if constexpr(I == 231) return (sqrt(104786045/(2.*pi))*y*(9*IntegerPow<8>(x) - 84*IntegerPow<6>(x)*IntegerPow<2>(y) + 126*IntegerPow<4>(x)*IntegerPow<4>(y) - 36*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*(-1 + 15*IntegerPow<2>(z)*(5 - 45*IntegerPow<2>(z) + 87*IntegerPow<4>(z))))/8192.;
    if constexpr(I == 232) return -(sqrt(44908305/pi)*y*(-IntegerPow<7>(x) + 7*IntegerPow<5>(x)*IntegerPow<2>(y) - 7*IntegerPow<3>(x)*IntegerPow<4>(y) + x*IntegerPow<6>(y))*z*(-7 + 175*IntegerPow<2>(z) - 945*IntegerPow<4>(z) + 1305*IntegerPow<6>(z)))/512.;
    if constexpr(I == 233) return -(sqrt(1952535/(2.*pi))*y*(-7*IntegerPow<6>(x) + 35*IntegerPow<4>(x)*IntegerPow<2>(y) - 21*IntegerPow<2>(x)*IntegerPow<4>(y) + IntegerPow<6>(y))*(7 + 23*IntegerPow<2>(z)*(-28 + 5*IntegerPow<2>(z)*(70 - 252*IntegerPow<2>(z) + 261*IntegerPow<4>(z)))))/8192.;
    if constexpr(I == 234) return (sqrt(21477885/pi)*y*(3*IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 3*x*IntegerPow<4>(y))*z*(21 + 23*IntegerPow<2>(z)*(-28 + 210*IntegerPow<2>(z) - 540*IntegerPow<4>(z) + 435*IntegerPow<6>(z))))/4096.;
    if constexpr(I == 235) return (3*sqrt(10023013/(2.*pi))*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(-1 + 5*IntegerPow<2>(z)*(21 + 23*IntegerPow<2>(z)*(-14 + 70*IntegerPow<2>(z) - 135*IntegerPow<4>(z) + 87*IntegerPow<6>(z)))))/8192.;
    if constexpr(I == 236) return (-3*sqrt(4555915/(2.*pi))*x*y*z*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(-11 + 385*IntegerPow<2>(z) + 23*IntegerPow<4>(z)*(-154 + 550*IntegerPow<2>(z) - 825*IntegerPow<4>(z) + 435*IntegerPow<6>(z))))/1024.;
    if constexpr(I == 237) return -(sqrt(719355/(2.*pi))*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*(11 + 19*IntegerPow<2>(z)*(-66 + 1155*IntegerPow<2>(z) + 23*IntegerPow<4>(z)*(-308 + 825*IntegerPow<2>(z) - 990*IntegerPow<4>(z) + 435*IntegerPow<6>(z)))))/8192.;
    if constexpr(I == 238) return (sqrt(55335/pi)*x*y*z*(429 + 19*IntegerPow<2>(z)*(-858 + 9009*IntegerPow<2>(z) + 23*IntegerPow<4>(z)*(-1716 + 5*IntegerPow<2>(z)*(715 - 702*IntegerPow<2>(z) + 261*IntegerPow<4>(z))))))/4096.;
    if constexpr(I == 239) return (sqrt(465/(2.*pi))*y*(-429 + 51051*IntegerPow<2>(z) + 323*IntegerPow<4>(z)*(-3003 + 21021*IntegerPow<2>(z) + 23*IntegerPow<4>(z)*(-3003 + 5005*IntegerPow<2>(z) - 4095*IntegerPow<4>(z) + 1305*IntegerPow<6>(z)))))/8192.;
    if constexpr(I == 240) return (sqrt(31/pi)*z*(-6435 + 255255*IntegerPow<2>(z) + 323*IntegerPow<4>(z)*(-9009 + 5*IntegerPow<2>(z)*(9009 + 23*IntegerPow<2>(z)*(-1001 + 1365*IntegerPow<2>(z) - 945*IntegerPow<4>(z) + 261*IntegerPow<6>(z))))))/4096.;
    if constexpr(I == 241) return (sqrt(465/(2.*pi))*x*(-429 + 51051*IntegerPow<2>(z) + 323*IntegerPow<4>(z)*(-3003 + 21021*IntegerPow<2>(z) + 23*IntegerPow<4>(z)*(-3003 + 5005*IntegerPow<2>(z) - 4095*IntegerPow<4>(z) + 1305*IntegerPow<6>(z)))))/8192.;
    if constexpr(I == 242) return -(sqrt(55335/pi)*z*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(429 + 19*IntegerPow<2>(z)*(-858 + 9009*IntegerPow<2>(z) + 23*IntegerPow<4>(z)*(-1716 + 5*IntegerPow<2>(z)*(715 - 702*IntegerPow<2>(z) + 261*IntegerPow<4>(z))))))/8192.;
    if constexpr(I == 243) return (sqrt(719355/(2.*pi))*(IntegerPow<3>(x) - 3*x*IntegerPow<2>(y))*(11 + 19*IntegerPow<2>(z)*(-66 + 1155*IntegerPow<2>(z) + 23*IntegerPow<4>(z)*(-308 + 825*IntegerPow<2>(z) - 990*IntegerPow<4>(z) + 435*IntegerPow<6>(z)))))/8192.;
    if constexpr(I == 244) return (3*sqrt(4555915/(2.*pi))*z*(8*IntegerPow<4>(y) + 8*IntegerPow<2>(y)*(-1 + IntegerPow<2>(z)) + IntegerPow<2>(-1 + IntegerPow<2>(z)))*(-11 + 385*IntegerPow<2>(z) + 23*IntegerPow<4>(z)*(-154 + 550*IntegerPow<2>(z) - 825*IntegerPow<4>(z) + 435*IntegerPow<6>(z))))/4096.;
    if constexpr(I == 245) return (3*sqrt(10023013/(2.*pi))*(IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 5*x*IntegerPow<4>(y))*(-1 + 5*IntegerPow<2>(z)*(21 + 23*IntegerPow<2>(z)*(-14 + 70*IntegerPow<2>(z) - 135*IntegerPow<4>(z) + 87*IntegerPow<6>(z)))))/8192.;
    if constexpr(I == 246) return -(sqrt(21477885/pi)*z*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(16*IntegerPow<4>(y) + 16*IntegerPow<2>(y)*(-1 + IntegerPow<2>(z)) + IntegerPow<2>(-1 + IntegerPow<2>(z)))*(21 + 23*IntegerPow<2>(z)*(-28 + 210*IntegerPow<2>(z) - 540*IntegerPow<4>(z) + 435*IntegerPow<6>(z))))/8192.;
    if constexpr(I == 247) return (sqrt(1952535/(2.*pi))*(IntegerPow<7>(x) - 21*IntegerPow<5>(x)*IntegerPow<2>(y) + 35*IntegerPow<3>(x)*IntegerPow<4>(y) - 7*x*IntegerPow<6>(y))*(7 + 23*IntegerPow<2>(z)*(-28 + 5*IntegerPow<2>(z)*(70 - 252*IntegerPow<2>(z) + 261*IntegerPow<4>(z)))))/8192.;
    if constexpr(I == 248) return (sqrt(44908305/pi)*(IntegerPow<8>(x) - 28*IntegerPow<6>(x)*IntegerPow<2>(y) + 70*IntegerPow<4>(x)*IntegerPow<4>(y) - 28*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*z*(-7 + 175*IntegerPow<2>(z) - 945*IntegerPow<4>(z) + 1305*IntegerPow<6>(z)))/4096.;
    if constexpr(I == 249) return (sqrt(104786045/(2.*pi))*(IntegerPow<9>(x) - 36*IntegerPow<7>(x)*IntegerPow<2>(y) + 126*IntegerPow<5>(x)*IntegerPow<4>(y) - 84*IntegerPow<3>(x)*IntegerPow<6>(y) + 9*x*IntegerPow<8>(y))*(-1 + 15*IntegerPow<2>(z)*(5 - 45*IntegerPow<2>(z) + 87*IntegerPow<4>(z))))/8192.;
    if constexpr(I == 250) return -(sqrt(314358135/pi)*z*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(5 - 90*IntegerPow<2>(z) + 261*IntegerPow<4>(z))*(256*IntegerPow<8>(y) + 512*IntegerPow<6>(y)*(-1 + IntegerPow<2>(z)) + 304*IntegerPow<4>(y)*IntegerPow<2>(-1 + IntegerPow<2>(z)) + 48*IntegerPow<2>(y)*IntegerPow<3>(-1 + IntegerPow<2>(z)) + IntegerPow<4>(-1 + IntegerPow<2>(z))))/8192.;
    if constexpr(I == 251) return (5*sqrt(4836279/(2.*pi))*(IntegerPow<11>(x) - 55*IntegerPow<9>(x)*IntegerPow<2>(y) + 330*IntegerPow<7>(x)*IntegerPow<4>(y) - 462*IntegerPow<5>(x)*IntegerPow<6>(y) + 165*IntegerPow<3>(x)*IntegerPow<8>(y) - 11*x*IntegerPow<10>(y))*(1 - 54*IntegerPow<2>(z) + 261*IntegerPow<4>(z)))/8192.;
    if constexpr(I == 252) return (15*sqrt(1612093/(2.*pi))*(IntegerPow<12>(x) - 66*IntegerPow<10>(x)*IntegerPow<2>(y) + 495*IntegerPow<8>(x)*IntegerPow<4>(y) - 924*IntegerPow<6>(x)*IntegerPow<6>(y) + 495*IntegerPow<4>(x)*IntegerPow<8>(y) - 66*IntegerPow<2>(x)*IntegerPow<10>(y) + IntegerPow<12>(y))*z*(-3 + 29*IntegerPow<2>(z)))/4096.;
    if constexpr(I == 253) return (15*sqrt(690897/(2.*pi))*(IntegerPow<13>(x) - 78*IntegerPow<11>(x)*IntegerPow<2>(y) + 715*IntegerPow<9>(x)*IntegerPow<4>(y) - 1716*IntegerPow<7>(x)*IntegerPow<6>(y) + 1287*IntegerPow<5>(x)*IntegerPow<8>(y) - 286*IntegerPow<3>(x)*IntegerPow<10>(y) + 13*x*IntegerPow<12>(y))*(-1 + 29*IntegerPow<2>(z)))/8192.;
    if constexpr(I == 254) return (15*sqrt(20036013/pi)*(IntegerPow<14>(x) - 91*IntegerPow<12>(x)*IntegerPow<2>(y) + 1001*IntegerPow<10>(x)*IntegerPow<4>(y) - 3003*IntegerPow<8>(x)*IntegerPow<6>(y) + 3003*IntegerPow<6>(x)*IntegerPow<8>(y) - 1001*IntegerPow<4>(x)*IntegerPow<10>(y) + 91*IntegerPow<2>(x)*IntegerPow<12>(y) - IntegerPow<14>(y))*z)/8192.;
    if constexpr(I == 255) return (3*sqrt(33393355/(2.*pi))*x*(IntegerPow<2>(x) - 3*IntegerPow<2>(y))*(IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + 5*IntegerPow<4>(y))*(IntegerPow<8>(x) - 92*IntegerPow<6>(x)*IntegerPow<2>(y) + 134*IntegerPow<4>(x)*IntegerPow<4>(y) - 28*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y)))/8192.;

    if constexpr(I > 255)
    {
        exit_on_error("Spherical Harmonic index i to big.");
        return 0;
    }
}
std::vector<double (*)(double, double, double)> SphericalHarmonicsXyz::Ylist =
{ &SphericalHarmonicsXyz::Y<  0>, &SphericalHarmonicsXyz::Y<  1>, &SphericalHarmonicsXyz::Y<  2>, &SphericalHarmonicsXyz::Y<  3>, &SphericalHarmonicsXyz::Y<  4>, &SphericalHarmonicsXyz::Y<  5>, &SphericalHarmonicsXyz::Y<  6>, &SphericalHarmonicsXyz::Y<  7>, &SphericalHarmonicsXyz::Y<  8>, &SphericalHarmonicsXyz::Y<  9>,
  &SphericalHarmonicsXyz::Y< 10>, &SphericalHarmonicsXyz::Y< 11>, &SphericalHarmonicsXyz::Y< 12>, &SphericalHarmonicsXyz::Y< 13>, &SphericalHarmonicsXyz::Y< 14>, &SphericalHarmonicsXyz::Y< 15>, &SphericalHarmonicsXyz::Y< 16>, &SphericalHarmonicsXyz::Y< 17>, &SphericalHarmonicsXyz::Y< 18>, &SphericalHarmonicsXyz::Y< 19>,
  &SphericalHarmonicsXyz::Y< 20>, &SphericalHarmonicsXyz::Y< 21>, &SphericalHarmonicsXyz::Y< 22>, &SphericalHarmonicsXyz::Y< 23>, &SphericalHarmonicsXyz::Y< 24>, &SphericalHarmonicsXyz::Y< 25>, &SphericalHarmonicsXyz::Y< 26>, &SphericalHarmonicsXyz::Y< 27>, &SphericalHarmonicsXyz::Y< 28>, &SphericalHarmonicsXyz::Y< 29>,
  &SphericalHarmonicsXyz::Y< 30>, &SphericalHarmonicsXyz::Y< 31>, &SphericalHarmonicsXyz::Y< 32>, &SphericalHarmonicsXyz::Y< 33>, &SphericalHarmonicsXyz::Y< 34>, &SphericalHarmonicsXyz::Y< 35>, &SphericalHarmonicsXyz::Y< 36>, &SphericalHarmonicsXyz::Y< 37>, &SphericalHarmonicsXyz::Y< 38>, &SphericalHarmonicsXyz::Y< 39>,
  &SphericalHarmonicsXyz::Y< 40>, &SphericalHarmonicsXyz::Y< 41>, &SphericalHarmonicsXyz::Y< 42>, &SphericalHarmonicsXyz::Y< 43>, &SphericalHarmonicsXyz::Y< 44>, &SphericalHarmonicsXyz::Y< 45>, &SphericalHarmonicsXyz::Y< 46>, &SphericalHarmonicsXyz::Y< 47>, &SphericalHarmonicsXyz::Y< 48>, &SphericalHarmonicsXyz::Y< 49>,
  &SphericalHarmonicsXyz::Y< 50>, &SphericalHarmonicsXyz::Y< 51>, &SphericalHarmonicsXyz::Y< 52>, &SphericalHarmonicsXyz::Y< 53>, &SphericalHarmonicsXyz::Y< 54>, &SphericalHarmonicsXyz::Y< 55>, &SphericalHarmonicsXyz::Y< 56>, &SphericalHarmonicsXyz::Y< 57>, &SphericalHarmonicsXyz::Y< 58>, &SphericalHarmonicsXyz::Y< 59>,
  &SphericalHarmonicsXyz::Y< 60>, &SphericalHarmonicsXyz::Y< 61>, &SphericalHarmonicsXyz::Y< 62>, &SphericalHarmonicsXyz::Y< 63>, &SphericalHarmonicsXyz::Y< 64>, &SphericalHarmonicsXyz::Y< 65>, &SphericalHarmonicsXyz::Y< 66>, &SphericalHarmonicsXyz::Y< 67>, &SphericalHarmonicsXyz::Y< 68>, &SphericalHarmonicsXyz::Y< 69>,
  &SphericalHarmonicsXyz::Y< 70>, &SphericalHarmonicsXyz::Y< 71>, &SphericalHarmonicsXyz::Y< 72>, &SphericalHarmonicsXyz::Y< 73>, &SphericalHarmonicsXyz::Y< 74>, &SphericalHarmonicsXyz::Y< 75>, &SphericalHarmonicsXyz::Y< 76>, &SphericalHarmonicsXyz::Y< 77>, &SphericalHarmonicsXyz::Y< 78>, &SphericalHarmonicsXyz::Y< 79>,
  &SphericalHarmonicsXyz::Y< 80>, &SphericalHarmonicsXyz::Y< 81>, &SphericalHarmonicsXyz::Y< 82>, &SphericalHarmonicsXyz::Y< 83>, &SphericalHarmonicsXyz::Y< 84>, &SphericalHarmonicsXyz::Y< 85>, &SphericalHarmonicsXyz::Y< 86>, &SphericalHarmonicsXyz::Y< 87>, &SphericalHarmonicsXyz::Y< 88>, &SphericalHarmonicsXyz::Y< 89>,
  &SphericalHarmonicsXyz::Y< 90>, &SphericalHarmonicsXyz::Y< 91>, &SphericalHarmonicsXyz::Y< 92>, &SphericalHarmonicsXyz::Y< 93>, &SphericalHarmonicsXyz::Y< 94>, &SphericalHarmonicsXyz::Y< 95>, &SphericalHarmonicsXyz::Y< 96>, &SphericalHarmonicsXyz::Y< 97>, &SphericalHarmonicsXyz::Y< 98>, &SphericalHarmonicsXyz::Y< 99>,
  &SphericalHarmonicsXyz::Y<100>, &SphericalHarmonicsXyz::Y<101>, &SphericalHarmonicsXyz::Y<102>, &SphericalHarmonicsXyz::Y<103>, &SphericalHarmonicsXyz::Y<104>, &SphericalHarmonicsXyz::Y<105>, &SphericalHarmonicsXyz::Y<106>, &SphericalHarmonicsXyz::Y<107>, &SphericalHarmonicsXyz::Y<108>, &SphericalHarmonicsXyz::Y<109>,
  &SphericalHarmonicsXyz::Y<110>, &SphericalHarmonicsXyz::Y<111>, &SphericalHarmonicsXyz::Y<112>, &SphericalHarmonicsXyz::Y<113>, &SphericalHarmonicsXyz::Y<114>, &SphericalHarmonicsXyz::Y<115>, &SphericalHarmonicsXyz::Y<116>, &SphericalHarmonicsXyz::Y<117>, &SphericalHarmonicsXyz::Y<118>, &SphericalHarmonicsXyz::Y<119>,
  &SphericalHarmonicsXyz::Y<120>, &SphericalHarmonicsXyz::Y<121>, &SphericalHarmonicsXyz::Y<122>, &SphericalHarmonicsXyz::Y<123>, &SphericalHarmonicsXyz::Y<124>, &SphericalHarmonicsXyz::Y<125>, &SphericalHarmonicsXyz::Y<126>, &SphericalHarmonicsXyz::Y<127>, &SphericalHarmonicsXyz::Y<128>, &SphericalHarmonicsXyz::Y<129>,
  &SphericalHarmonicsXyz::Y<130>, &SphericalHarmonicsXyz::Y<131>, &SphericalHarmonicsXyz::Y<132>, &SphericalHarmonicsXyz::Y<133>, &SphericalHarmonicsXyz::Y<134>, &SphericalHarmonicsXyz::Y<135>, &SphericalHarmonicsXyz::Y<136>, &SphericalHarmonicsXyz::Y<137>, &SphericalHarmonicsXyz::Y<138>, &SphericalHarmonicsXyz::Y<139>,
  &SphericalHarmonicsXyz::Y<140>, &SphericalHarmonicsXyz::Y<141>, &SphericalHarmonicsXyz::Y<142>, &SphericalHarmonicsXyz::Y<143>, &SphericalHarmonicsXyz::Y<144>, &SphericalHarmonicsXyz::Y<145>, &SphericalHarmonicsXyz::Y<146>, &SphericalHarmonicsXyz::Y<147>, &SphericalHarmonicsXyz::Y<148>, &SphericalHarmonicsXyz::Y<149>,
  &SphericalHarmonicsXyz::Y<150>, &SphericalHarmonicsXyz::Y<151>, &SphericalHarmonicsXyz::Y<152>, &SphericalHarmonicsXyz::Y<153>, &SphericalHarmonicsXyz::Y<154>, &SphericalHarmonicsXyz::Y<155>, &SphericalHarmonicsXyz::Y<156>, &SphericalHarmonicsXyz::Y<157>, &SphericalHarmonicsXyz::Y<158>, &SphericalHarmonicsXyz::Y<159>,
  &SphericalHarmonicsXyz::Y<160>, &SphericalHarmonicsXyz::Y<161>, &SphericalHarmonicsXyz::Y<162>, &SphericalHarmonicsXyz::Y<163>, &SphericalHarmonicsXyz::Y<164>, &SphericalHarmonicsXyz::Y<165>, &SphericalHarmonicsXyz::Y<166>, &SphericalHarmonicsXyz::Y<167>, &SphericalHarmonicsXyz::Y<168>, &SphericalHarmonicsXyz::Y<169>,
  &SphericalHarmonicsXyz::Y<170>, &SphericalHarmonicsXyz::Y<171>, &SphericalHarmonicsXyz::Y<172>, &SphericalHarmonicsXyz::Y<173>, &SphericalHarmonicsXyz::Y<174>, &SphericalHarmonicsXyz::Y<175>, &SphericalHarmonicsXyz::Y<176>, &SphericalHarmonicsXyz::Y<177>, &SphericalHarmonicsXyz::Y<178>, &SphericalHarmonicsXyz::Y<179>,
  &SphericalHarmonicsXyz::Y<180>, &SphericalHarmonicsXyz::Y<181>, &SphericalHarmonicsXyz::Y<182>, &SphericalHarmonicsXyz::Y<183>, &SphericalHarmonicsXyz::Y<184>, &SphericalHarmonicsXyz::Y<185>, &SphericalHarmonicsXyz::Y<186>, &SphericalHarmonicsXyz::Y<187>, &SphericalHarmonicsXyz::Y<188>, &SphericalHarmonicsXyz::Y<189>,
  &SphericalHarmonicsXyz::Y<190>, &SphericalHarmonicsXyz::Y<191>, &SphericalHarmonicsXyz::Y<192>, &SphericalHarmonicsXyz::Y<193>, &SphericalHarmonicsXyz::Y<194>, &SphericalHarmonicsXyz::Y<195>, &SphericalHarmonicsXyz::Y<196>, &SphericalHarmonicsXyz::Y<197>, &SphericalHarmonicsXyz::Y<198>, &SphericalHarmonicsXyz::Y<199>,
  &SphericalHarmonicsXyz::Y<200>, &SphericalHarmonicsXyz::Y<201>, &SphericalHarmonicsXyz::Y<202>, &SphericalHarmonicsXyz::Y<203>, &SphericalHarmonicsXyz::Y<204>, &SphericalHarmonicsXyz::Y<205>, &SphericalHarmonicsXyz::Y<206>, &SphericalHarmonicsXyz::Y<207>, &SphericalHarmonicsXyz::Y<208>, &SphericalHarmonicsXyz::Y<209>,
  &SphericalHarmonicsXyz::Y<210>, &SphericalHarmonicsXyz::Y<211>, &SphericalHarmonicsXyz::Y<212>, &SphericalHarmonicsXyz::Y<213>, &SphericalHarmonicsXyz::Y<214>, &SphericalHarmonicsXyz::Y<215>, &SphericalHarmonicsXyz::Y<216>, &SphericalHarmonicsXyz::Y<217>, &SphericalHarmonicsXyz::Y<218>, &SphericalHarmonicsXyz::Y<219>,
  &SphericalHarmonicsXyz::Y<220>, &SphericalHarmonicsXyz::Y<221>, &SphericalHarmonicsXyz::Y<222>, &SphericalHarmonicsXyz::Y<223>, &SphericalHarmonicsXyz::Y<224>, &SphericalHarmonicsXyz::Y<225>, &SphericalHarmonicsXyz::Y<226>, &SphericalHarmonicsXyz::Y<227>, &SphericalHarmonicsXyz::Y<228>, &SphericalHarmonicsXyz::Y<229>,
  &SphericalHarmonicsXyz::Y<230>, &SphericalHarmonicsXyz::Y<231>, &SphericalHarmonicsXyz::Y<232>, &SphericalHarmonicsXyz::Y<233>, &SphericalHarmonicsXyz::Y<234>, &SphericalHarmonicsXyz::Y<235>, &SphericalHarmonicsXyz::Y<236>, &SphericalHarmonicsXyz::Y<237>, &SphericalHarmonicsXyz::Y<238>, &SphericalHarmonicsXyz::Y<239>,
  &SphericalHarmonicsXyz::Y<240>, &SphericalHarmonicsXyz::Y<241>, &SphericalHarmonicsXyz::Y<242>, &SphericalHarmonicsXyz::Y<243>, &SphericalHarmonicsXyz::Y<244>, &SphericalHarmonicsXyz::Y<245>, &SphericalHarmonicsXyz::Y<246>, &SphericalHarmonicsXyz::Y<247>, &SphericalHarmonicsXyz::Y<248>, &SphericalHarmonicsXyz::Y<249>,
  &SphericalHarmonicsXyz::Y<250>, &SphericalHarmonicsXyz::Y<251>, &SphericalHarmonicsXyz::Y<252>, &SphericalHarmonicsXyz::Y<253>, &SphericalHarmonicsXyz::Y<254>, &SphericalHarmonicsXyz::Y<255> };



// Y_i(x,y,z), i=l(l+1)+m.
double SphericalHarmonicsXyz::Y(int i, double x, double y, double z)
{
    // i = l*(l+1) + m, with l € [0,N], m € [-l,l]
    // l can be reconstructed by taking the sqrt of i and rounding down:    int l = floor(sqrt(i));
    // m can then simply be obtained by rearranged the above equation:      int m = i - l * (l + 1);
    return Ylist[i](x,y,z);
}
double SphericalHarmonicsXyz::Y(int i, const Tensor3& dir)
{
    // i = l*(l+1) + m, with l € [0,N], m € [-l,l]
    // l can be reconstructed by taking the sqrt of i and rounding down:    int l = floor(sqrt(i));
    // m can then simply be obtained by rearranged the above equation:      int m = i - l * (l + 1);
    return Ylist[i](dir[1],dir[2],dir[3]);
}
// Y_lm(x,y,z).
double SphericalHarmonicsXyz::Y(int l, int m, double x, double y, double z)
{
    // l € [0,N]
    // m € [-l,l]
    // theta € [0, pi]
    //   phi € [0,2pi]
    int i = l * (l + 1) + m;
    return Ylist[i](x,y,z);
}
double SphericalHarmonicsXyz::Y(int l, int m, const Tensor3& dir)
{
    // l € [0,N]
    // m € [-l,l]
    // theta € [0, pi]
    //   phi € [0,2pi]
    int i = l * (l + 1) + m;
    return Ylist[i](dir[1],dir[2],dir[3]);
}



std::vector<double> SphericalHarmonicsXyz::GetCoefficients(const Stencil& stencil, const double* data)
{
    std::vector<double> coefficients(stencil.nCoefficients);
    for(size_t d=0; d<stencil.nDir; d++)
    {
        double x = stencil.Cx(d);
        double y = stencil.Cy(d);
        double z = stencil.Cz(d);
        double c = data[d] * stencil.W(d);

        for(size_t i=0; i<stencil.nCoefficients; i++)
            coefficients[i] += c * Y(i,x,y,z);
    }
    return coefficients;
}
void SphericalHarmonicsXyz::GetCoefficients(const Stencil& stencil, const double* data, double* coefficients)
{
    for(size_t i=0; i<stencil.nCoefficients; i++)
        coefficients[i] = 0;
    for(size_t d=0; d<stencil.nDir; d++)
    {
        double x = stencil.Cx(d);
        double y = stencil.Cy(d);
        double z = stencil.Cz(d);
        double c = data[d] * stencil.W(d);

        for(size_t i=0; i<stencil.nCoefficients; i++)
            coefficients[i] += c * Y(i,x,y,z);
    }
}
double SphericalHarmonicsXyz::GetValue(double x, double y, double z, const std::vector<double>& coefficients, size_t nCoefficients)
{
    double result = 0;
    for(size_t i=0; i<nCoefficients; i++)
        result += coefficients[i] * Y(i,x,y,z);

    return result;
}
double SphericalHarmonicsXyz::GetValue(double x, double y, double z, double* coefficients, size_t nCoefficients)
{
    double result = 0;
    for(size_t i=0; i<nCoefficients; i++)
        result += coefficients[i] * Y(i,x,y,z);

    return result;
}
double SphericalHarmonicsXyz::GetValue(const Tensor3& direction, const std::vector<double>& coefficients, size_t nCoefficients)
{
    double result = 0;
    for(size_t i=0; i<nCoefficients; i++)
        result += coefficients[i] * Y(i,direction[1],direction[2],direction[3]);

    return result;
}
double SphericalHarmonicsXyz::GetValue(const Tensor3& direction, double* coefficients, size_t nCoefficients)
{
    double result = 0;
    for(size_t i=0; i<nCoefficients; i++)
        result += coefficients[i] * Y(i,direction[1],direction[2],direction[3]);

    return result;
}
// ---------------------------------------------------------------------------------------------