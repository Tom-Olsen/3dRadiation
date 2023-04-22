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
    if constexpr(I == 256) return Y<16,-16>(theta, phi);
    if constexpr(I == 257) return Y<16,-15>(theta, phi);
    if constexpr(I == 258) return Y<16,-14>(theta, phi);
    if constexpr(I == 259) return Y<16,-13>(theta, phi);
    if constexpr(I == 260) return Y<16,-12>(theta, phi);
    if constexpr(I == 261) return Y<16,-11>(theta, phi);
    if constexpr(I == 262) return Y<16,-10>(theta, phi);
    if constexpr(I == 263) return Y<16, -9>(theta, phi);
    if constexpr(I == 264) return Y<16, -8>(theta, phi);
    if constexpr(I == 265) return Y<16, -7>(theta, phi);
    if constexpr(I == 266) return Y<16, -6>(theta, phi);
    if constexpr(I == 267) return Y<16, -5>(theta, phi);
    if constexpr(I == 268) return Y<16, -4>(theta, phi);
    if constexpr(I == 269) return Y<16, -3>(theta, phi);
    if constexpr(I == 270) return Y<16, -2>(theta, phi);
    if constexpr(I == 271) return Y<16, -1>(theta, phi);
    if constexpr(I == 272) return Y<16,  0>(theta, phi);
    if constexpr(I == 273) return Y<16,  1>(theta, phi);
    if constexpr(I == 274) return Y<16,  2>(theta, phi);
    if constexpr(I == 275) return Y<16,  3>(theta, phi);
    if constexpr(I == 276) return Y<16,  4>(theta, phi);
    if constexpr(I == 277) return Y<16,  5>(theta, phi);
    if constexpr(I == 278) return Y<16,  6>(theta, phi);
    if constexpr(I == 279) return Y<16,  7>(theta, phi);
    if constexpr(I == 280) return Y<16,  8>(theta, phi);
    if constexpr(I == 281) return Y<16,  9>(theta, phi);
    if constexpr(I == 282) return Y<16, 10>(theta, phi);
    if constexpr(I == 283) return Y<16, 11>(theta, phi);
    if constexpr(I == 284) return Y<16, 12>(theta, phi);
    if constexpr(I == 285) return Y<16, 13>(theta, phi);
    if constexpr(I == 286) return Y<16, 14>(theta, phi);
    if constexpr(I == 287) return Y<16, 15>(theta, phi);
    if constexpr(I == 288) return Y<16, 16>(theta, phi);
    if constexpr(I == 289) return Y<17,-17>(theta, phi);
    if constexpr(I == 290) return Y<17,-16>(theta, phi);
    if constexpr(I == 291) return Y<17,-15>(theta, phi);
    if constexpr(I == 292) return Y<17,-14>(theta, phi);
    if constexpr(I == 293) return Y<17,-13>(theta, phi);
    if constexpr(I == 294) return Y<17,-12>(theta, phi);
    if constexpr(I == 295) return Y<17,-11>(theta, phi);
    if constexpr(I == 296) return Y<17,-10>(theta, phi);
    if constexpr(I == 297) return Y<17, -9>(theta, phi);
    if constexpr(I == 298) return Y<17, -8>(theta, phi);
    if constexpr(I == 299) return Y<17, -7>(theta, phi);
    if constexpr(I == 300) return Y<17, -6>(theta, phi);
    if constexpr(I == 301) return Y<17, -5>(theta, phi);
    if constexpr(I == 302) return Y<17, -4>(theta, phi);
    if constexpr(I == 303) return Y<17, -3>(theta, phi);
    if constexpr(I == 304) return Y<17, -2>(theta, phi);
    if constexpr(I == 305) return Y<17, -1>(theta, phi);
    if constexpr(I == 306) return Y<17,  0>(theta, phi);
    if constexpr(I == 307) return Y<17,  1>(theta, phi);
    if constexpr(I == 308) return Y<17,  2>(theta, phi);
    if constexpr(I == 309) return Y<17,  3>(theta, phi);
    if constexpr(I == 310) return Y<17,  4>(theta, phi);
    if constexpr(I == 311) return Y<17,  5>(theta, phi);
    if constexpr(I == 312) return Y<17,  6>(theta, phi);
    if constexpr(I == 313) return Y<17,  7>(theta, phi);
    if constexpr(I == 314) return Y<17,  8>(theta, phi);
    if constexpr(I == 315) return Y<17,  9>(theta, phi);
    if constexpr(I == 316) return Y<17, 10>(theta, phi);
    if constexpr(I == 317) return Y<17, 11>(theta, phi);
    if constexpr(I == 318) return Y<17, 12>(theta, phi);
    if constexpr(I == 319) return Y<17, 13>(theta, phi);
    if constexpr(I == 320) return Y<17, 14>(theta, phi);
    if constexpr(I == 321) return Y<17, 15>(theta, phi);
    if constexpr(I == 322) return Y<17, 16>(theta, phi);
    if constexpr(I == 323) return Y<17, 17>(theta, phi);
    if constexpr(I == 324) return Y<18,-18>(theta, phi);
    if constexpr(I == 325) return Y<18,-17>(theta, phi);
    if constexpr(I == 326) return Y<18,-16>(theta, phi);
    if constexpr(I == 327) return Y<18,-15>(theta, phi);
    if constexpr(I == 328) return Y<18,-14>(theta, phi);
    if constexpr(I == 329) return Y<18,-13>(theta, phi);
    if constexpr(I == 330) return Y<18,-12>(theta, phi);
    if constexpr(I == 331) return Y<18,-11>(theta, phi);
    if constexpr(I == 332) return Y<18,-10>(theta, phi);
    if constexpr(I == 333) return Y<18, -9>(theta, phi);
    if constexpr(I == 334) return Y<18, -8>(theta, phi);
    if constexpr(I == 335) return Y<18, -7>(theta, phi);
    if constexpr(I == 336) return Y<18, -6>(theta, phi);
    if constexpr(I == 337) return Y<18, -5>(theta, phi);
    if constexpr(I == 338) return Y<18, -4>(theta, phi);
    if constexpr(I == 339) return Y<18, -3>(theta, phi);
    if constexpr(I == 340) return Y<18, -2>(theta, phi);
    if constexpr(I == 341) return Y<18, -1>(theta, phi);
    if constexpr(I == 342) return Y<18,  0>(theta, phi);
    if constexpr(I == 343) return Y<18,  1>(theta, phi);
    if constexpr(I == 344) return Y<18,  2>(theta, phi);
    if constexpr(I == 345) return Y<18,  3>(theta, phi);
    if constexpr(I == 346) return Y<18,  4>(theta, phi);
    if constexpr(I == 347) return Y<18,  5>(theta, phi);
    if constexpr(I == 348) return Y<18,  6>(theta, phi);
    if constexpr(I == 349) return Y<18,  7>(theta, phi);
    if constexpr(I == 350) return Y<18,  8>(theta, phi);
    if constexpr(I == 351) return Y<18,  9>(theta, phi);
    if constexpr(I == 352) return Y<18, 10>(theta, phi);
    if constexpr(I == 353) return Y<18, 11>(theta, phi);
    if constexpr(I == 354) return Y<18, 12>(theta, phi);
    if constexpr(I == 355) return Y<18, 13>(theta, phi);
    if constexpr(I == 356) return Y<18, 14>(theta, phi);
    if constexpr(I == 357) return Y<18, 15>(theta, phi);
    if constexpr(I == 358) return Y<18, 16>(theta, phi);
    if constexpr(I == 359) return Y<18, 17>(theta, phi);
    if constexpr(I == 360) return Y<18, 18>(theta, phi);
    if constexpr(I == 361) return Y<19,-19>(theta, phi);
    if constexpr(I == 362) return Y<19,-18>(theta, phi);
    if constexpr(I == 363) return Y<19,-17>(theta, phi);
    if constexpr(I == 364) return Y<19,-16>(theta, phi);
    if constexpr(I == 365) return Y<19,-15>(theta, phi);
    if constexpr(I == 366) return Y<19,-14>(theta, phi);
    if constexpr(I == 367) return Y<19,-13>(theta, phi);
    if constexpr(I == 368) return Y<19,-12>(theta, phi);
    if constexpr(I == 369) return Y<19,-11>(theta, phi);
    if constexpr(I == 370) return Y<19,-10>(theta, phi);
    if constexpr(I == 371) return Y<19, -9>(theta, phi);
    if constexpr(I == 372) return Y<19, -8>(theta, phi);
    if constexpr(I == 373) return Y<19, -7>(theta, phi);
    if constexpr(I == 374) return Y<19, -6>(theta, phi);
    if constexpr(I == 375) return Y<19, -5>(theta, phi);
    if constexpr(I == 376) return Y<19, -4>(theta, phi);
    if constexpr(I == 377) return Y<19, -3>(theta, phi);
    if constexpr(I == 378) return Y<19, -2>(theta, phi);
    if constexpr(I == 379) return Y<19, -1>(theta, phi);
    if constexpr(I == 380) return Y<19,  0>(theta, phi);
    if constexpr(I == 381) return Y<19,  1>(theta, phi);
    if constexpr(I == 382) return Y<19,  2>(theta, phi);
    if constexpr(I == 383) return Y<19,  3>(theta, phi);
    if constexpr(I == 384) return Y<19,  4>(theta, phi);
    if constexpr(I == 385) return Y<19,  5>(theta, phi);
    if constexpr(I == 386) return Y<19,  6>(theta, phi);
    if constexpr(I == 387) return Y<19,  7>(theta, phi);
    if constexpr(I == 388) return Y<19,  8>(theta, phi);
    if constexpr(I == 389) return Y<19,  9>(theta, phi);
    if constexpr(I == 390) return Y<19, 10>(theta, phi);
    if constexpr(I == 391) return Y<19, 11>(theta, phi);
    if constexpr(I == 392) return Y<19, 12>(theta, phi);
    if constexpr(I == 393) return Y<19, 13>(theta, phi);
    if constexpr(I == 394) return Y<19, 14>(theta, phi);
    if constexpr(I == 395) return Y<19, 15>(theta, phi);
    if constexpr(I == 396) return Y<19, 16>(theta, phi);
    if constexpr(I == 397) return Y<19, 17>(theta, phi);
    if constexpr(I == 398) return Y<19, 18>(theta, phi);
    if constexpr(I == 399) return Y<19, 19>(theta, phi);
    if constexpr(I == 400) return Y<20,-20>(theta, phi);
    if constexpr(I == 401) return Y<20,-19>(theta, phi);
    if constexpr(I == 402) return Y<20,-18>(theta, phi);
    if constexpr(I == 403) return Y<20,-17>(theta, phi);
    if constexpr(I == 404) return Y<20,-16>(theta, phi);
    if constexpr(I == 405) return Y<20,-15>(theta, phi);
    if constexpr(I == 406) return Y<20,-14>(theta, phi);
    if constexpr(I == 407) return Y<20,-13>(theta, phi);
    if constexpr(I == 408) return Y<20,-12>(theta, phi);
    if constexpr(I == 409) return Y<20,-11>(theta, phi);
    if constexpr(I == 410) return Y<20,-10>(theta, phi);
    if constexpr(I == 411) return Y<20, -9>(theta, phi);
    if constexpr(I == 412) return Y<20, -8>(theta, phi);
    if constexpr(I == 413) return Y<20, -7>(theta, phi);
    if constexpr(I == 414) return Y<20, -6>(theta, phi);
    if constexpr(I == 415) return Y<20, -5>(theta, phi);
    if constexpr(I == 416) return Y<20, -4>(theta, phi);
    if constexpr(I == 417) return Y<20, -3>(theta, phi);
    if constexpr(I == 418) return Y<20, -2>(theta, phi);
    if constexpr(I == 419) return Y<20, -1>(theta, phi);
    if constexpr(I == 420) return Y<20,  0>(theta, phi);
    if constexpr(I == 421) return Y<20,  1>(theta, phi);
    if constexpr(I == 422) return Y<20,  2>(theta, phi);
    if constexpr(I == 423) return Y<20,  3>(theta, phi);
    if constexpr(I == 424) return Y<20,  4>(theta, phi);
    if constexpr(I == 425) return Y<20,  5>(theta, phi);
    if constexpr(I == 426) return Y<20,  6>(theta, phi);
    if constexpr(I == 427) return Y<20,  7>(theta, phi);
    if constexpr(I == 428) return Y<20,  8>(theta, phi);
    if constexpr(I == 429) return Y<20,  9>(theta, phi);
    if constexpr(I == 430) return Y<20, 10>(theta, phi);
    if constexpr(I == 431) return Y<20, 11>(theta, phi);
    if constexpr(I == 432) return Y<20, 12>(theta, phi);
    if constexpr(I == 433) return Y<20, 13>(theta, phi);
    if constexpr(I == 434) return Y<20, 14>(theta, phi);
    if constexpr(I == 435) return Y<20, 15>(theta, phi);
    if constexpr(I == 436) return Y<20, 16>(theta, phi);
    if constexpr(I == 437) return Y<20, 17>(theta, phi);
    if constexpr(I == 438) return Y<20, 18>(theta, phi);
    if constexpr(I == 439) return Y<20, 19>(theta, phi);
    if constexpr(I == 440) return Y<20, 20>(theta, phi);
    if constexpr(I == 441) return Y<21,-21>(theta, phi);
    if constexpr(I == 442) return Y<21,-20>(theta, phi);
    if constexpr(I == 443) return Y<21,-19>(theta, phi);
    if constexpr(I == 444) return Y<21,-18>(theta, phi);
    if constexpr(I == 445) return Y<21,-17>(theta, phi);
    if constexpr(I == 446) return Y<21,-16>(theta, phi);
    if constexpr(I == 447) return Y<21,-15>(theta, phi);
    if constexpr(I == 448) return Y<21,-14>(theta, phi);
    if constexpr(I == 449) return Y<21,-13>(theta, phi);
    if constexpr(I == 450) return Y<21,-12>(theta, phi);
    if constexpr(I == 451) return Y<21,-11>(theta, phi);
    if constexpr(I == 452) return Y<21,-10>(theta, phi);
    if constexpr(I == 453) return Y<21, -9>(theta, phi);
    if constexpr(I == 454) return Y<21, -8>(theta, phi);
    if constexpr(I == 455) return Y<21, -7>(theta, phi);
    if constexpr(I == 456) return Y<21, -6>(theta, phi);
    if constexpr(I == 457) return Y<21, -5>(theta, phi);
    if constexpr(I == 458) return Y<21, -4>(theta, phi);
    if constexpr(I == 459) return Y<21, -3>(theta, phi);
    if constexpr(I == 460) return Y<21, -2>(theta, phi);
    if constexpr(I == 461) return Y<21, -1>(theta, phi);
    if constexpr(I == 462) return Y<21,  0>(theta, phi);
    if constexpr(I == 463) return Y<21,  1>(theta, phi);
    if constexpr(I == 464) return Y<21,  2>(theta, phi);
    if constexpr(I == 465) return Y<21,  3>(theta, phi);
    if constexpr(I == 466) return Y<21,  4>(theta, phi);
    if constexpr(I == 467) return Y<21,  5>(theta, phi);
    if constexpr(I == 468) return Y<21,  6>(theta, phi);
    if constexpr(I == 469) return Y<21,  7>(theta, phi);
    if constexpr(I == 470) return Y<21,  8>(theta, phi);
    if constexpr(I == 471) return Y<21,  9>(theta, phi);
    if constexpr(I == 472) return Y<21, 10>(theta, phi);
    if constexpr(I == 473) return Y<21, 11>(theta, phi);
    if constexpr(I == 474) return Y<21, 12>(theta, phi);
    if constexpr(I == 475) return Y<21, 13>(theta, phi);
    if constexpr(I == 476) return Y<21, 14>(theta, phi);
    if constexpr(I == 477) return Y<21, 15>(theta, phi);
    if constexpr(I == 478) return Y<21, 16>(theta, phi);
    if constexpr(I == 479) return Y<21, 17>(theta, phi);
    if constexpr(I == 480) return Y<21, 18>(theta, phi);
    if constexpr(I == 481) return Y<21, 19>(theta, phi);
    if constexpr(I == 482) return Y<21, 20>(theta, phi);
    if constexpr(I == 483) return Y<21, 21>(theta, phi);
    if constexpr(I == 484) return Y<22,-22>(theta, phi);
    if constexpr(I == 485) return Y<22,-21>(theta, phi);
    if constexpr(I == 486) return Y<22,-20>(theta, phi);
    if constexpr(I == 487) return Y<22,-19>(theta, phi);
    if constexpr(I == 488) return Y<22,-18>(theta, phi);
    if constexpr(I == 489) return Y<22,-17>(theta, phi);
    if constexpr(I == 490) return Y<22,-16>(theta, phi);
    if constexpr(I == 491) return Y<22,-15>(theta, phi);
    if constexpr(I == 492) return Y<22,-14>(theta, phi);
    if constexpr(I == 493) return Y<22,-13>(theta, phi);
    if constexpr(I == 494) return Y<22,-12>(theta, phi);
    if constexpr(I == 495) return Y<22,-11>(theta, phi);
    if constexpr(I == 496) return Y<22,-10>(theta, phi);
    if constexpr(I == 497) return Y<22, -9>(theta, phi);
    if constexpr(I == 498) return Y<22, -8>(theta, phi);
    if constexpr(I == 499) return Y<22, -7>(theta, phi);
    if constexpr(I == 500) return Y<22, -6>(theta, phi);
    if constexpr(I == 501) return Y<22, -5>(theta, phi);
    if constexpr(I == 502) return Y<22, -4>(theta, phi);
    if constexpr(I == 503) return Y<22, -3>(theta, phi);
    if constexpr(I == 504) return Y<22, -2>(theta, phi);
    if constexpr(I == 505) return Y<22, -1>(theta, phi);
    if constexpr(I == 506) return Y<22,  0>(theta, phi);
    if constexpr(I == 507) return Y<22,  1>(theta, phi);
    if constexpr(I == 508) return Y<22,  2>(theta, phi);
    if constexpr(I == 509) return Y<22,  3>(theta, phi);
    if constexpr(I == 510) return Y<22,  4>(theta, phi);
    if constexpr(I == 511) return Y<22,  5>(theta, phi);
    if constexpr(I == 512) return Y<22,  6>(theta, phi);
    if constexpr(I == 513) return Y<22,  7>(theta, phi);
    if constexpr(I == 514) return Y<22,  8>(theta, phi);
    if constexpr(I == 515) return Y<22,  9>(theta, phi);
    if constexpr(I == 516) return Y<22, 10>(theta, phi);
    if constexpr(I == 517) return Y<22, 11>(theta, phi);
    if constexpr(I == 518) return Y<22, 12>(theta, phi);
    if constexpr(I == 519) return Y<22, 13>(theta, phi);
    if constexpr(I == 520) return Y<22, 14>(theta, phi);
    if constexpr(I == 521) return Y<22, 15>(theta, phi);
    if constexpr(I == 522) return Y<22, 16>(theta, phi);
    if constexpr(I == 523) return Y<22, 17>(theta, phi);
    if constexpr(I == 524) return Y<22, 18>(theta, phi);
    if constexpr(I == 525) return Y<22, 19>(theta, phi);
    if constexpr(I == 526) return Y<22, 20>(theta, phi);
    if constexpr(I == 527) return Y<22, 21>(theta, phi);
    if constexpr(I == 528) return Y<22, 22>(theta, phi);
    if constexpr(I == 529) return Y<23,-23>(theta, phi);
    if constexpr(I == 530) return Y<23,-22>(theta, phi);
    if constexpr(I == 531) return Y<23,-21>(theta, phi);
    if constexpr(I == 532) return Y<23,-20>(theta, phi);
    if constexpr(I == 533) return Y<23,-19>(theta, phi);
    if constexpr(I == 534) return Y<23,-18>(theta, phi);
    if constexpr(I == 535) return Y<23,-17>(theta, phi);
    if constexpr(I == 536) return Y<23,-16>(theta, phi);
    if constexpr(I == 537) return Y<23,-15>(theta, phi);
    if constexpr(I == 538) return Y<23,-14>(theta, phi);
    if constexpr(I == 539) return Y<23,-13>(theta, phi);
    if constexpr(I == 540) return Y<23,-12>(theta, phi);
    if constexpr(I == 541) return Y<23,-11>(theta, phi);
    if constexpr(I == 542) return Y<23,-10>(theta, phi);
    if constexpr(I == 543) return Y<23, -9>(theta, phi);
    if constexpr(I == 544) return Y<23, -8>(theta, phi);
    if constexpr(I == 545) return Y<23, -7>(theta, phi);
    if constexpr(I == 546) return Y<23, -6>(theta, phi);
    if constexpr(I == 547) return Y<23, -5>(theta, phi);
    if constexpr(I == 548) return Y<23, -4>(theta, phi);
    if constexpr(I == 549) return Y<23, -3>(theta, phi);
    if constexpr(I == 550) return Y<23, -2>(theta, phi);
    if constexpr(I == 551) return Y<23, -1>(theta, phi);
    if constexpr(I == 552) return Y<23,  0>(theta, phi);
    if constexpr(I == 553) return Y<23,  1>(theta, phi);
    if constexpr(I == 554) return Y<23,  2>(theta, phi);
    if constexpr(I == 555) return Y<23,  3>(theta, phi);
    if constexpr(I == 556) return Y<23,  4>(theta, phi);
    if constexpr(I == 557) return Y<23,  5>(theta, phi);
    if constexpr(I == 558) return Y<23,  6>(theta, phi);
    if constexpr(I == 559) return Y<23,  7>(theta, phi);
    if constexpr(I == 560) return Y<23,  8>(theta, phi);
    if constexpr(I == 561) return Y<23,  9>(theta, phi);
    if constexpr(I == 562) return Y<23, 10>(theta, phi);
    if constexpr(I == 563) return Y<23, 11>(theta, phi);
    if constexpr(I == 564) return Y<23, 12>(theta, phi);
    if constexpr(I == 565) return Y<23, 13>(theta, phi);
    if constexpr(I == 566) return Y<23, 14>(theta, phi);
    if constexpr(I == 567) return Y<23, 15>(theta, phi);
    if constexpr(I == 568) return Y<23, 16>(theta, phi);
    if constexpr(I == 569) return Y<23, 17>(theta, phi);
    if constexpr(I == 570) return Y<23, 18>(theta, phi);
    if constexpr(I == 571) return Y<23, 19>(theta, phi);
    if constexpr(I == 572) return Y<23, 20>(theta, phi);
    if constexpr(I == 573) return Y<23, 21>(theta, phi);
    if constexpr(I == 574) return Y<23, 22>(theta, phi);
    if constexpr(I == 575) return Y<23, 23>(theta, phi);
    if constexpr(I > 575)
    {
        ExitOnError("Spherical Harmonic index i to big.");
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
  &SphericalHarmonicsThPh::Y<250>, &SphericalHarmonicsThPh::Y<251>, &SphericalHarmonicsThPh::Y<252>, &SphericalHarmonicsThPh::Y<253>, &SphericalHarmonicsThPh::Y<254>, &SphericalHarmonicsThPh::Y<255>, &SphericalHarmonicsThPh::Y<256>, &SphericalHarmonicsThPh::Y<257>, &SphericalHarmonicsThPh::Y<258>, &SphericalHarmonicsThPh::Y<259>,
  &SphericalHarmonicsThPh::Y<260>, &SphericalHarmonicsThPh::Y<261>, &SphericalHarmonicsThPh::Y<262>, &SphericalHarmonicsThPh::Y<263>, &SphericalHarmonicsThPh::Y<264>, &SphericalHarmonicsThPh::Y<265>, &SphericalHarmonicsThPh::Y<266>, &SphericalHarmonicsThPh::Y<267>, &SphericalHarmonicsThPh::Y<268>, &SphericalHarmonicsThPh::Y<269>,
  &SphericalHarmonicsThPh::Y<270>, &SphericalHarmonicsThPh::Y<271>, &SphericalHarmonicsThPh::Y<272>, &SphericalHarmonicsThPh::Y<273>, &SphericalHarmonicsThPh::Y<274>, &SphericalHarmonicsThPh::Y<275>, &SphericalHarmonicsThPh::Y<276>, &SphericalHarmonicsThPh::Y<277>, &SphericalHarmonicsThPh::Y<278>, &SphericalHarmonicsThPh::Y<279>,
  &SphericalHarmonicsThPh::Y<280>, &SphericalHarmonicsThPh::Y<281>, &SphericalHarmonicsThPh::Y<282>, &SphericalHarmonicsThPh::Y<283>, &SphericalHarmonicsThPh::Y<284>, &SphericalHarmonicsThPh::Y<285>, &SphericalHarmonicsThPh::Y<286>, &SphericalHarmonicsThPh::Y<287>, &SphericalHarmonicsThPh::Y<288>, &SphericalHarmonicsThPh::Y<289>,
  &SphericalHarmonicsThPh::Y<290>, &SphericalHarmonicsThPh::Y<291>, &SphericalHarmonicsThPh::Y<292>, &SphericalHarmonicsThPh::Y<293>, &SphericalHarmonicsThPh::Y<294>, &SphericalHarmonicsThPh::Y<295>, &SphericalHarmonicsThPh::Y<296>, &SphericalHarmonicsThPh::Y<297>, &SphericalHarmonicsThPh::Y<298>, &SphericalHarmonicsThPh::Y<299>,
  &SphericalHarmonicsThPh::Y<300>, &SphericalHarmonicsThPh::Y<301>, &SphericalHarmonicsThPh::Y<302>, &SphericalHarmonicsThPh::Y<303>, &SphericalHarmonicsThPh::Y<304>, &SphericalHarmonicsThPh::Y<305>, &SphericalHarmonicsThPh::Y<306>, &SphericalHarmonicsThPh::Y<307>, &SphericalHarmonicsThPh::Y<308>, &SphericalHarmonicsThPh::Y<309>,
  &SphericalHarmonicsThPh::Y<310>, &SphericalHarmonicsThPh::Y<311>, &SphericalHarmonicsThPh::Y<312>, &SphericalHarmonicsThPh::Y<313>, &SphericalHarmonicsThPh::Y<314>, &SphericalHarmonicsThPh::Y<315>, &SphericalHarmonicsThPh::Y<316>, &SphericalHarmonicsThPh::Y<317>, &SphericalHarmonicsThPh::Y<318>, &SphericalHarmonicsThPh::Y<319>,
  &SphericalHarmonicsThPh::Y<320>, &SphericalHarmonicsThPh::Y<321>, &SphericalHarmonicsThPh::Y<322>, &SphericalHarmonicsThPh::Y<323>, &SphericalHarmonicsThPh::Y<324>, &SphericalHarmonicsThPh::Y<325>, &SphericalHarmonicsThPh::Y<326>, &SphericalHarmonicsThPh::Y<327>, &SphericalHarmonicsThPh::Y<328>, &SphericalHarmonicsThPh::Y<329>,
  &SphericalHarmonicsThPh::Y<330>, &SphericalHarmonicsThPh::Y<331>, &SphericalHarmonicsThPh::Y<332>, &SphericalHarmonicsThPh::Y<333>, &SphericalHarmonicsThPh::Y<334>, &SphericalHarmonicsThPh::Y<335>, &SphericalHarmonicsThPh::Y<336>, &SphericalHarmonicsThPh::Y<337>, &SphericalHarmonicsThPh::Y<338>, &SphericalHarmonicsThPh::Y<339>,
  &SphericalHarmonicsThPh::Y<340>, &SphericalHarmonicsThPh::Y<341>, &SphericalHarmonicsThPh::Y<342>, &SphericalHarmonicsThPh::Y<343>, &SphericalHarmonicsThPh::Y<344>, &SphericalHarmonicsThPh::Y<345>, &SphericalHarmonicsThPh::Y<346>, &SphericalHarmonicsThPh::Y<347>, &SphericalHarmonicsThPh::Y<348>, &SphericalHarmonicsThPh::Y<349>,
  &SphericalHarmonicsThPh::Y<350>, &SphericalHarmonicsThPh::Y<351>, &SphericalHarmonicsThPh::Y<352>, &SphericalHarmonicsThPh::Y<353>, &SphericalHarmonicsThPh::Y<354>, &SphericalHarmonicsThPh::Y<355>, &SphericalHarmonicsThPh::Y<356>, &SphericalHarmonicsThPh::Y<357>, &SphericalHarmonicsThPh::Y<358>, &SphericalHarmonicsThPh::Y<359>,
  &SphericalHarmonicsThPh::Y<360>, &SphericalHarmonicsThPh::Y<361>, &SphericalHarmonicsThPh::Y<362>, &SphericalHarmonicsThPh::Y<363>, &SphericalHarmonicsThPh::Y<364>, &SphericalHarmonicsThPh::Y<365>, &SphericalHarmonicsThPh::Y<366>, &SphericalHarmonicsThPh::Y<367>, &SphericalHarmonicsThPh::Y<368>, &SphericalHarmonicsThPh::Y<369>,
  &SphericalHarmonicsThPh::Y<370>, &SphericalHarmonicsThPh::Y<371>, &SphericalHarmonicsThPh::Y<372>, &SphericalHarmonicsThPh::Y<373>, &SphericalHarmonicsThPh::Y<374>, &SphericalHarmonicsThPh::Y<375>, &SphericalHarmonicsThPh::Y<376>, &SphericalHarmonicsThPh::Y<377>, &SphericalHarmonicsThPh::Y<378>, &SphericalHarmonicsThPh::Y<379>,
  &SphericalHarmonicsThPh::Y<380>, &SphericalHarmonicsThPh::Y<381>, &SphericalHarmonicsThPh::Y<382>, &SphericalHarmonicsThPh::Y<383>, &SphericalHarmonicsThPh::Y<384>, &SphericalHarmonicsThPh::Y<385>, &SphericalHarmonicsThPh::Y<386>, &SphericalHarmonicsThPh::Y<387>, &SphericalHarmonicsThPh::Y<388>, &SphericalHarmonicsThPh::Y<389>,
  &SphericalHarmonicsThPh::Y<390>, &SphericalHarmonicsThPh::Y<391>, &SphericalHarmonicsThPh::Y<392>, &SphericalHarmonicsThPh::Y<393>, &SphericalHarmonicsThPh::Y<394>, &SphericalHarmonicsThPh::Y<395>, &SphericalHarmonicsThPh::Y<396>, &SphericalHarmonicsThPh::Y<397>, &SphericalHarmonicsThPh::Y<398>, &SphericalHarmonicsThPh::Y<399>,
  &SphericalHarmonicsThPh::Y<400>, &SphericalHarmonicsThPh::Y<401>, &SphericalHarmonicsThPh::Y<402>, &SphericalHarmonicsThPh::Y<403>, &SphericalHarmonicsThPh::Y<404>, &SphericalHarmonicsThPh::Y<405>, &SphericalHarmonicsThPh::Y<406>, &SphericalHarmonicsThPh::Y<407>, &SphericalHarmonicsThPh::Y<408>, &SphericalHarmonicsThPh::Y<409>,
  &SphericalHarmonicsThPh::Y<410>, &SphericalHarmonicsThPh::Y<411>, &SphericalHarmonicsThPh::Y<412>, &SphericalHarmonicsThPh::Y<413>, &SphericalHarmonicsThPh::Y<414>, &SphericalHarmonicsThPh::Y<415>, &SphericalHarmonicsThPh::Y<416>, &SphericalHarmonicsThPh::Y<417>, &SphericalHarmonicsThPh::Y<418>, &SphericalHarmonicsThPh::Y<419>,
  &SphericalHarmonicsThPh::Y<420>, &SphericalHarmonicsThPh::Y<421>, &SphericalHarmonicsThPh::Y<422>, &SphericalHarmonicsThPh::Y<423>, &SphericalHarmonicsThPh::Y<424>, &SphericalHarmonicsThPh::Y<425>, &SphericalHarmonicsThPh::Y<426>, &SphericalHarmonicsThPh::Y<427>, &SphericalHarmonicsThPh::Y<428>, &SphericalHarmonicsThPh::Y<429>,
  &SphericalHarmonicsThPh::Y<430>, &SphericalHarmonicsThPh::Y<431>, &SphericalHarmonicsThPh::Y<432>, &SphericalHarmonicsThPh::Y<433>, &SphericalHarmonicsThPh::Y<434>, &SphericalHarmonicsThPh::Y<435>, &SphericalHarmonicsThPh::Y<436>, &SphericalHarmonicsThPh::Y<437>, &SphericalHarmonicsThPh::Y<438>, &SphericalHarmonicsThPh::Y<439>,
  &SphericalHarmonicsThPh::Y<440>, &SphericalHarmonicsThPh::Y<441>, &SphericalHarmonicsThPh::Y<442>, &SphericalHarmonicsThPh::Y<443>, &SphericalHarmonicsThPh::Y<444>, &SphericalHarmonicsThPh::Y<445>, &SphericalHarmonicsThPh::Y<446>, &SphericalHarmonicsThPh::Y<447>, &SphericalHarmonicsThPh::Y<448>, &SphericalHarmonicsThPh::Y<449>,
  &SphericalHarmonicsThPh::Y<450>, &SphericalHarmonicsThPh::Y<451>, &SphericalHarmonicsThPh::Y<452>, &SphericalHarmonicsThPh::Y<453>, &SphericalHarmonicsThPh::Y<454>, &SphericalHarmonicsThPh::Y<455>, &SphericalHarmonicsThPh::Y<456>, &SphericalHarmonicsThPh::Y<457>, &SphericalHarmonicsThPh::Y<458>, &SphericalHarmonicsThPh::Y<459>,
  &SphericalHarmonicsThPh::Y<460>, &SphericalHarmonicsThPh::Y<461>, &SphericalHarmonicsThPh::Y<462>, &SphericalHarmonicsThPh::Y<463>, &SphericalHarmonicsThPh::Y<464>, &SphericalHarmonicsThPh::Y<465>, &SphericalHarmonicsThPh::Y<466>, &SphericalHarmonicsThPh::Y<467>, &SphericalHarmonicsThPh::Y<468>, &SphericalHarmonicsThPh::Y<469>,
  &SphericalHarmonicsThPh::Y<470>, &SphericalHarmonicsThPh::Y<471>, &SphericalHarmonicsThPh::Y<472>, &SphericalHarmonicsThPh::Y<473>, &SphericalHarmonicsThPh::Y<474>, &SphericalHarmonicsThPh::Y<475>, &SphericalHarmonicsThPh::Y<476>, &SphericalHarmonicsThPh::Y<477>, &SphericalHarmonicsThPh::Y<478>, &SphericalHarmonicsThPh::Y<479>,
  &SphericalHarmonicsThPh::Y<480>, &SphericalHarmonicsThPh::Y<481>, &SphericalHarmonicsThPh::Y<482>, &SphericalHarmonicsThPh::Y<483>, &SphericalHarmonicsThPh::Y<484>, &SphericalHarmonicsThPh::Y<485>, &SphericalHarmonicsThPh::Y<486>, &SphericalHarmonicsThPh::Y<487>, &SphericalHarmonicsThPh::Y<488>, &SphericalHarmonicsThPh::Y<489>,
  &SphericalHarmonicsThPh::Y<490>, &SphericalHarmonicsThPh::Y<491>, &SphericalHarmonicsThPh::Y<492>, &SphericalHarmonicsThPh::Y<493>, &SphericalHarmonicsThPh::Y<494>, &SphericalHarmonicsThPh::Y<495>, &SphericalHarmonicsThPh::Y<496>, &SphericalHarmonicsThPh::Y<497>, &SphericalHarmonicsThPh::Y<498>, &SphericalHarmonicsThPh::Y<499>,
  &SphericalHarmonicsThPh::Y<500>, &SphericalHarmonicsThPh::Y<501>, &SphericalHarmonicsThPh::Y<502>, &SphericalHarmonicsThPh::Y<503>, &SphericalHarmonicsThPh::Y<504>, &SphericalHarmonicsThPh::Y<505>, &SphericalHarmonicsThPh::Y<506>, &SphericalHarmonicsThPh::Y<507>, &SphericalHarmonicsThPh::Y<508>, &SphericalHarmonicsThPh::Y<509>,
  &SphericalHarmonicsThPh::Y<510>, &SphericalHarmonicsThPh::Y<511>, &SphericalHarmonicsThPh::Y<512>, &SphericalHarmonicsThPh::Y<513>, &SphericalHarmonicsThPh::Y<514>, &SphericalHarmonicsThPh::Y<515>, &SphericalHarmonicsThPh::Y<516>, &SphericalHarmonicsThPh::Y<517>, &SphericalHarmonicsThPh::Y<518>, &SphericalHarmonicsThPh::Y<519>,
  &SphericalHarmonicsThPh::Y<520>, &SphericalHarmonicsThPh::Y<521>, &SphericalHarmonicsThPh::Y<522>, &SphericalHarmonicsThPh::Y<523>, &SphericalHarmonicsThPh::Y<524>, &SphericalHarmonicsThPh::Y<525>, &SphericalHarmonicsThPh::Y<526>, &SphericalHarmonicsThPh::Y<527>, &SphericalHarmonicsThPh::Y<528>, &SphericalHarmonicsThPh::Y<529>,
  &SphericalHarmonicsThPh::Y<530>, &SphericalHarmonicsThPh::Y<531>, &SphericalHarmonicsThPh::Y<532>, &SphericalHarmonicsThPh::Y<533>, &SphericalHarmonicsThPh::Y<534>, &SphericalHarmonicsThPh::Y<535>, &SphericalHarmonicsThPh::Y<536>, &SphericalHarmonicsThPh::Y<537>, &SphericalHarmonicsThPh::Y<538>, &SphericalHarmonicsThPh::Y<539>,
  &SphericalHarmonicsThPh::Y<540>, &SphericalHarmonicsThPh::Y<541>, &SphericalHarmonicsThPh::Y<542>, &SphericalHarmonicsThPh::Y<543>, &SphericalHarmonicsThPh::Y<544>, &SphericalHarmonicsThPh::Y<545>, &SphericalHarmonicsThPh::Y<546>, &SphericalHarmonicsThPh::Y<547>, &SphericalHarmonicsThPh::Y<548>, &SphericalHarmonicsThPh::Y<549>,
  &SphericalHarmonicsThPh::Y<550>, &SphericalHarmonicsThPh::Y<551>, &SphericalHarmonicsThPh::Y<552>, &SphericalHarmonicsThPh::Y<553>, &SphericalHarmonicsThPh::Y<554>, &SphericalHarmonicsThPh::Y<555>, &SphericalHarmonicsThPh::Y<556>, &SphericalHarmonicsThPh::Y<557>, &SphericalHarmonicsThPh::Y<558>, &SphericalHarmonicsThPh::Y<559>,
  &SphericalHarmonicsThPh::Y<560>, &SphericalHarmonicsThPh::Y<561>, &SphericalHarmonicsThPh::Y<562>, &SphericalHarmonicsThPh::Y<563>, &SphericalHarmonicsThPh::Y<564>, &SphericalHarmonicsThPh::Y<565>, &SphericalHarmonicsThPh::Y<566>, &SphericalHarmonicsThPh::Y<567>, &SphericalHarmonicsThPh::Y<568>, &SphericalHarmonicsThPh::Y<569>,
  &SphericalHarmonicsThPh::Y<570>, &SphericalHarmonicsThPh::Y<571>, &SphericalHarmonicsThPh::Y<572>, &SphericalHarmonicsThPh::Y<573>, &SphericalHarmonicsThPh::Y<574>, &SphericalHarmonicsThPh::Y<575> };



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
    if constexpr(I == 256) return (3*sqrt(1101980715/pi)*x*y*(IntegerPow<14>(x) - 35*IntegerPow<12>(x)*IntegerPow<2>(y) + 273*IntegerPow<10>(x)*IntegerPow<4>(y) - 715*IntegerPow<8>(x)*IntegerPow<6>(y) + 715*IntegerPow<6>(x)*IntegerPow<8>(y) - 273*IntegerPow<4>(x)*IntegerPow<10>(y) + 35*IntegerPow<2>(x)*IntegerPow<12>(y) - IntegerPow<14>(y)))/4096.;
    if constexpr(I == 257) return (-3*sqrt(1101980715/(2.*pi))*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(IntegerPow<8>(x) - 28*IntegerPow<6>(x)*IntegerPow<2>(y) + 134*IntegerPow<4>(x)*IntegerPow<4>(y) - 92*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*z)/8192.;
    if constexpr(I == 258) return (3*sqrt(35547765/pi)*x*y*(IntegerPow<6>(x) - 21*IntegerPow<4>(x)*IntegerPow<2>(y) + 35*IntegerPow<2>(x)*IntegerPow<4>(y) - 7*IntegerPow<6>(y))*(7*IntegerPow<6>(x) - 35*IntegerPow<4>(x)*IntegerPow<2>(y) + 21*IntegerPow<2>(x)*IntegerPow<4>(y) - IntegerPow<6>(y))*(-1 + 31*IntegerPow<2>(z)))/8192.;
    if constexpr(I == 259) return (15*sqrt(7109553/(2.*pi))*y*(13*IntegerPow<12>(x) - 286*IntegerPow<10>(x)*IntegerPow<2>(y) + 1287*IntegerPow<8>(x)*IntegerPow<4>(y) - 1716*IntegerPow<6>(x)*IntegerPow<6>(y) + 715*IntegerPow<4>(x)*IntegerPow<8>(y) - 78*IntegerPow<2>(x)*IntegerPow<10>(y) + IntegerPow<12>(y))*z*(-3 + 31*IntegerPow<2>(z)))/8192.;
    if constexpr(I == 260) return (15*sqrt(245157/(2.*pi))*x*(x - y)*y*(x + y)*(IntegerPow<2>(x) - 3*IntegerPow<2>(y))*(3*IntegerPow<2>(x) - IntegerPow<2>(y))*(IntegerPow<2>(x) - 4*x*y + IntegerPow<2>(y))*(IntegerPow<2>(x) + 4*x*y + IntegerPow<2>(y))*(3 - 174*IntegerPow<2>(z) + 899*IntegerPow<4>(z)))/4096.;
    if constexpr(I == 261) return (-3*sqrt(8580495/(2.*pi))*y*(-11*IntegerPow<10>(x) + 165*IntegerPow<8>(x)*IntegerPow<2>(y) - 462*IntegerPow<6>(x)*IntegerPow<4>(y) + 330*IntegerPow<4>(x)*IntegerPow<6>(y) - 55*IntegerPow<2>(x)*IntegerPow<8>(y) + IntegerPow<10>(y))*z*(15 - 290*IntegerPow<2>(z) + 899*IntegerPow<4>(z)))/8192.;
    if constexpr(I == 262) return (sqrt(8580495/pi)*x*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + 5*IntegerPow<4>(y))*(-5 + 9*IntegerPow<2>(z)*(45 - 435*IntegerPow<2>(z) + 899*IntegerPow<4>(z))))/8192.;
    if constexpr(I == 263) return (sqrt(15935205/(2.*pi))*y*(9*IntegerPow<8>(x) - 84*IntegerPow<6>(x)*IntegerPow<2>(y) + 126*IntegerPow<4>(x)*IntegerPow<4>(y) - 36*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*z*(-35 + 945*IntegerPow<2>(z) - 5481*IntegerPow<4>(z) + 8091*IntegerPow<6>(z)))/8192.;
    if constexpr(I == 264) return -(sqrt(15935205/pi)*y*(-IntegerPow<7>(x) + 7*IntegerPow<5>(x)*IntegerPow<2>(y) - 7*IntegerPow<3>(x)*IntegerPow<4>(y) + x*IntegerPow<6>(y))*(7 - 700*IntegerPow<2>(z) + 45*IntegerPow<4>(z)*(210 - 812*IntegerPow<2>(z) + 899*IntegerPow<4>(z))))/4096.;
    if constexpr(I == 265) return (-3*sqrt(5311735/(2.*pi))*y*(-7*IntegerPow<6>(x) + 35*IntegerPow<4>(x)*IntegerPow<2>(y) - 21*IntegerPow<2>(x)*IntegerPow<4>(y) + IntegerPow<6>(y))*z*(21 - 700*IntegerPow<2>(z) + 5670*IntegerPow<4>(z) - 15660*IntegerPow<6>(z) + 13485*IntegerPow<8>(z)))/8192.;
    if constexpr(I == 266) return (3*sqrt(46189/pi)*y*(3*IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 3*x*IntegerPow<4>(y))*(-21 + 115*IntegerPow<2>(z)*(21 - 350*IntegerPow<2>(z) + 1890*IntegerPow<4>(z) - 3915*IntegerPow<6>(z) + 2697*IntegerPow<8>(z))))/8192.;
    if constexpr(I == 267) return (3*sqrt(46189/(2.*pi))*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*z*(-231 + 115*IntegerPow<2>(z)*(77 - 770*IntegerPow<2>(z) + 2970*IntegerPow<4>(z) - 4785*IntegerPow<6>(z) + 2697*IntegerPow<8>(z))))/8192.;
    if constexpr(I == 268) return (-3*sqrt(323323/(2.*pi))*x*y*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(11 - 1386*IntegerPow<2>(z) + 115*IntegerPow<4>(z)*(231 - 1540*IntegerPow<2>(z) + 4455*IntegerPow<4>(z) - 5742*IntegerPow<6>(z) + 2697*IntegerPow<8>(z))))/4096.;
    if constexpr(I == 269) return (-3*sqrt(124355/(2.*pi))*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*z*(143 - 6006*IntegerPow<2>(z) + 23*IntegerPow<4>(z)*(3003 + 5*IntegerPow<2>(z)*(-2860 + 6435*IntegerPow<2>(z) - 6786*IntegerPow<4>(z) + 2697*IntegerPow<6>(z)))))/8192.;
    if constexpr(I == 270) return (3*sqrt(935/pi)*x*y*(-143 + 19*IntegerPow<2>(z)*(1001 - 21021*IntegerPow<2>(z) + 23*IntegerPow<4>(z)*(7007 + 5*IntegerPow<2>(z)*(-5005 + 9009*IntegerPow<2>(z) - 7917*IntegerPow<4>(z) + 2697*IntegerPow<6>(z))))))/8192.;
    if constexpr(I == 271) return (sqrt(561/(2.*pi))*y*z*(-6435 + 19*IntegerPow<2>(z)*(15015 - 189189*IntegerPow<2>(z) + 115*IntegerPow<4>(z)*(9009 - 25025*IntegerPow<2>(z) + 36855*IntegerPow<4>(z) - 27405*IntegerPow<6>(z) + 8091*IntegerPow<8>(z)))))/8192.;
    if constexpr(I == 272) return (sqrt(33/pi)*(6435 - 875160*IntegerPow<2>(z) + 323*IntegerPow<4>(z)*(60060 - 504504*IntegerPow<2>(z) + 115*IntegerPow<4>(z)*(18018 - 40040*IntegerPow<2>(z) + 49140*IntegerPow<4>(z) - 31320*IntegerPow<6>(z) + 8091*IntegerPow<8>(z)))))/65536.;
    if constexpr(I == 273) return (sqrt(561/(2.*pi))*x*z*(-6435 + 19*IntegerPow<2>(z)*(15015 - 189189*IntegerPow<2>(z) + 115*IntegerPow<4>(z)*(9009 - 25025*IntegerPow<2>(z) + 36855*IntegerPow<4>(z) - 27405*IntegerPow<6>(z) + 8091*IntegerPow<8>(z)))))/8192.;
    if constexpr(I == 274) return (3*sqrt(935/pi)*(x - y)*(x + y)*(-143 + 19*IntegerPow<2>(z)*(1001 - 21021*IntegerPow<2>(z) + 23*IntegerPow<4>(z)*(7007 + 5*IntegerPow<2>(z)*(-5005 + 9009*IntegerPow<2>(z) - 7917*IntegerPow<4>(z) + 2697*IntegerPow<6>(z))))))/16384.;
    if constexpr(I == 275) return (3*sqrt(124355/(2.*pi))*(IntegerPow<3>(x) - 3*x*IntegerPow<2>(y))*z*(143 - 6006*IntegerPow<2>(z) + 23*IntegerPow<4>(z)*(3003 + 5*IntegerPow<2>(z)*(-2860 + 6435*IntegerPow<2>(z) - 6786*IntegerPow<4>(z) + 2697*IntegerPow<6>(z)))))/8192.;
    if constexpr(I == 276) return (3*sqrt(323323/(2.*pi))*(8*IntegerPow<4>(y) + 8*IntegerPow<2>(y)*(-1 + IntegerPow<2>(z)) + IntegerPow<2>(-1 + IntegerPow<2>(z)))*(11 - 1386*IntegerPow<2>(z) + 115*IntegerPow<4>(z)*(231 - 1540*IntegerPow<2>(z) + 4455*IntegerPow<4>(z) - 5742*IntegerPow<6>(z) + 2697*IntegerPow<8>(z))))/16384.;
    if constexpr(I == 277) return (3*sqrt(46189/(2.*pi))*(IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 5*x*IntegerPow<4>(y))*z*(-231 + 115*IntegerPow<2>(z)*(77 - 770*IntegerPow<2>(z) + 2970*IntegerPow<4>(z) - 4785*IntegerPow<6>(z) + 2697*IntegerPow<8>(z))))/8192.;
    if constexpr(I == 278) return (-3*sqrt(46189/pi)*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(16*IntegerPow<4>(y) + 16*IntegerPow<2>(y)*(-1 + IntegerPow<2>(z)) + IntegerPow<2>(-1 + IntegerPow<2>(z)))*(-21 + 115*IntegerPow<2>(z)*(21 - 350*IntegerPow<2>(z) + 1890*IntegerPow<4>(z) - 3915*IntegerPow<6>(z) + 2697*IntegerPow<8>(z))))/16384.;
    if constexpr(I == 279) return (3*sqrt(5311735/(2.*pi))*(IntegerPow<7>(x) - 21*IntegerPow<5>(x)*IntegerPow<2>(y) + 35*IntegerPow<3>(x)*IntegerPow<4>(y) - 7*x*IntegerPow<6>(y))*z*(21 - 700*IntegerPow<2>(z) + 5670*IntegerPow<4>(z) - 15660*IntegerPow<6>(z) + 13485*IntegerPow<8>(z)))/8192.;
    if constexpr(I == 280) return (sqrt(15935205/pi)*(128*IntegerPow<8>(y) + 256*IntegerPow<6>(y)*(-1 + IntegerPow<2>(z)) + 160*IntegerPow<4>(y)*IntegerPow<2>(-1 + IntegerPow<2>(z)) + 32*IntegerPow<2>(y)*IntegerPow<3>(-1 + IntegerPow<2>(z)) + IntegerPow<4>(-1 + IntegerPow<2>(z)))*(7 - 700*IntegerPow<2>(z) + 45*IntegerPow<4>(z)*(210 - 812*IntegerPow<2>(z) + 899*IntegerPow<4>(z))))/32768.;
    if constexpr(I == 281) return (sqrt(15935205/(2.*pi))*(IntegerPow<9>(x) - 36*IntegerPow<7>(x)*IntegerPow<2>(y) + 126*IntegerPow<5>(x)*IntegerPow<4>(y) - 84*IntegerPow<3>(x)*IntegerPow<6>(y) + 9*x*IntegerPow<8>(y))*z*(-35 + 945*IntegerPow<2>(z) - 5481*IntegerPow<4>(z) + 8091*IntegerPow<6>(z)))/8192.;
    if constexpr(I == 282) return (sqrt(8580495/pi)*(IntegerPow<10>(x) - 45*IntegerPow<8>(x)*IntegerPow<2>(y) + 210*IntegerPow<6>(x)*IntegerPow<4>(y) - 210*IntegerPow<4>(x)*IntegerPow<6>(y) + 45*IntegerPow<2>(x)*IntegerPow<8>(y) - IntegerPow<10>(y))*(-5 + 9*IntegerPow<2>(z)*(45 - 435*IntegerPow<2>(z) + 899*IntegerPow<4>(z))))/16384.;
    if constexpr(I == 283) return (3*sqrt(8580495/(2.*pi))*(IntegerPow<11>(x) - 55*IntegerPow<9>(x)*IntegerPow<2>(y) + 330*IntegerPow<7>(x)*IntegerPow<4>(y) - 462*IntegerPow<5>(x)*IntegerPow<6>(y) + 165*IntegerPow<3>(x)*IntegerPow<8>(y) - 11*x*IntegerPow<10>(y))*z*(15 - 290*IntegerPow<2>(z) + 899*IntegerPow<4>(z)))/8192.;
    if constexpr(I == 284) return (15*sqrt(245157/(2.*pi))*(3 - 174*IntegerPow<2>(z) + 899*IntegerPow<4>(z))*(2048*IntegerPow<12>(y) + 6144*IntegerPow<10>(y)*(-1 + IntegerPow<2>(z)) + 6912*IntegerPow<8>(y)*IntegerPow<2>(-1 + IntegerPow<2>(z)) + 3584*IntegerPow<6>(y)*IntegerPow<3>(-1 + IntegerPow<2>(z)) + 840*IntegerPow<4>(y)*IntegerPow<4>(-1 + IntegerPow<2>(z)) + 72*IntegerPow<2>(y)*IntegerPow<5>(-1 + IntegerPow<2>(z)) + IntegerPow<6>(-1 + IntegerPow<2>(z))))/16384.;
    if constexpr(I == 285) return (15*sqrt(7109553/(2.*pi))*(IntegerPow<13>(x) - 78*IntegerPow<11>(x)*IntegerPow<2>(y) + 715*IntegerPow<9>(x)*IntegerPow<4>(y) - 1716*IntegerPow<7>(x)*IntegerPow<6>(y) + 1287*IntegerPow<5>(x)*IntegerPow<8>(y) - 286*IntegerPow<3>(x)*IntegerPow<10>(y) + 13*x*IntegerPow<12>(y))*z*(-3 + 31*IntegerPow<2>(z)))/8192.;
    if constexpr(I == 286) return (3*sqrt(35547765/pi)*(IntegerPow<14>(x) - 91*IntegerPow<12>(x)*IntegerPow<2>(y) + 1001*IntegerPow<10>(x)*IntegerPow<4>(y) - 3003*IntegerPow<8>(x)*IntegerPow<6>(y) + 3003*IntegerPow<6>(x)*IntegerPow<8>(y) - 1001*IntegerPow<4>(x)*IntegerPow<10>(y) + 91*IntegerPow<2>(x)*IntegerPow<12>(y) - IntegerPow<14>(y))*(-1 + 31*IntegerPow<2>(z)))/16384.;
    if constexpr(I == 287) return (3*sqrt(1101980715/(2.*pi))*x*(IntegerPow<2>(x) - 3*IntegerPow<2>(y))*(IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + 5*IntegerPow<4>(y))*(IntegerPow<8>(x) - 92*IntegerPow<6>(x)*IntegerPow<2>(y) + 134*IntegerPow<4>(x)*IntegerPow<4>(y) - 28*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*z)/8192.;
    if constexpr(I == 288) return (3*sqrt(1101980715/pi)*(IntegerPow<16>(x) - 120*IntegerPow<14>(x)*IntegerPow<2>(y) + 1820*IntegerPow<12>(x)*IntegerPow<4>(y) - 8008*IntegerPow<10>(x)*IntegerPow<6>(y) + 12870*IntegerPow<8>(x)*IntegerPow<8>(y) - 8008*IntegerPow<6>(x)*IntegerPow<10>(y) + 1820*IntegerPow<4>(x)*IntegerPow<12>(y) - 120*IntegerPow<2>(x)*IntegerPow<14>(y) + IntegerPow<16>(y)))/65536.;
    if constexpr(I == 289) return (15*sqrt(90751353/(2.*pi))*y*(17*IntegerPow<16>(x) - 680*IntegerPow<14>(x)*IntegerPow<2>(y) + 6188*IntegerPow<12>(x)*IntegerPow<4>(y) - 19448*IntegerPow<10>(x)*IntegerPow<6>(y) + 24310*IntegerPow<8>(x)*IntegerPow<8>(y) - 12376*IntegerPow<6>(x)*IntegerPow<10>(y) + 2380*IntegerPow<4>(x)*IntegerPow<12>(y) - 136*IntegerPow<2>(x)*IntegerPow<14>(y) + IntegerPow<16>(y)))/65536.;
    if constexpr(I == 290) return (15*sqrt(1542773001/pi)*x*y*(IntegerPow<14>(x) - 35*IntegerPow<12>(x)*IntegerPow<2>(y) + 273*IntegerPow<10>(x)*IntegerPow<4>(y) - 715*IntegerPow<8>(x)*IntegerPow<6>(y) + 715*IntegerPow<6>(x)*IntegerPow<8>(y) - 273*IntegerPow<4>(x)*IntegerPow<10>(y) + 35*IntegerPow<2>(x)*IntegerPow<12>(y) - IntegerPow<14>(y))*z)/4096.;
    if constexpr(I == 291) return (-15*sqrt(46750697/(2.*pi))*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(IntegerPow<8>(x) - 28*IntegerPow<6>(x)*IntegerPow<2>(y) + 134*IntegerPow<4>(x)*IntegerPow<4>(y) - 92*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*(-1 + 33*IntegerPow<2>(z)))/65536.;
    if constexpr(I == 292) return (15*sqrt(140252091/pi)*x*y*(IntegerPow<6>(x) - 21*IntegerPow<4>(x)*IntegerPow<2>(y) + 35*IntegerPow<2>(x)*IntegerPow<4>(y) - 7*IntegerPow<6>(y))*(7*IntegerPow<6>(x) - 35*IntegerPow<4>(x)*IntegerPow<2>(y) + 21*IntegerPow<2>(x)*IntegerPow<4>(y) - IntegerPow<6>(y))*z*(-1 + 11*IntegerPow<2>(z)))/8192.;
    if constexpr(I == 293) return (15*sqrt(4524261/pi)*y*(13*IntegerPow<12>(x) - 286*IntegerPow<10>(x)*IntegerPow<2>(y) + 1287*IntegerPow<8>(x)*IntegerPow<4>(y) - 1716*IntegerPow<6>(x)*IntegerPow<6>(y) + 715*IntegerPow<4>(x)*IntegerPow<8>(y) - 78*IntegerPow<2>(x)*IntegerPow<10>(y) + IntegerPow<12>(y))*(1 - 62*IntegerPow<2>(z) + 341*IntegerPow<4>(z)))/32768.;
    if constexpr(I == 294) return (15*sqrt(1508087/(2.*pi))*x*(x - y)*y*(x + y)*(IntegerPow<2>(x) - 3*IntegerPow<2>(y))*(3*IntegerPow<2>(x) - IntegerPow<2>(y))*(IntegerPow<2>(x) - 4*x*y + IntegerPow<2>(y))*(IntegerPow<2>(x) + 4*x*y + IntegerPow<2>(y))*z*(15 - 310*IntegerPow<2>(z) + 1023*IntegerPow<4>(z)))/4096.;
    if constexpr(I == 295) return (-15*sqrt(156009/pi)*y*(-11*IntegerPow<10>(x) + 165*IntegerPow<8>(x)*IntegerPow<2>(y) - 462*IntegerPow<6>(x)*IntegerPow<4>(y) + 330*IntegerPow<4>(x)*IntegerPow<6>(y) - 55*IntegerPow<2>(x)*IntegerPow<8>(y) + IntegerPow<10>(y))*(-5 + 435*IntegerPow<2>(z) - 4495*IntegerPow<4>(z) + 9889*IntegerPow<6>(z)))/32768.;
    if constexpr(I == 296) return (15*sqrt(156009/pi)*x*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + 5*IntegerPow<4>(y))*z*(-35 + 29*IntegerPow<2>(z)*(35 - 217*IntegerPow<2>(z) + 341*IntegerPow<4>(z))))/8192.;
    if constexpr(I == 297) return (5*sqrt(52003/(2.*pi))*y*(9*IntegerPow<8>(x) - 84*IntegerPow<6>(x)*IntegerPow<2>(y) + 126*IntegerPow<4>(x)*IntegerPow<4>(y) - 36*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*(35 + 9*IntegerPow<2>(z)*(-420 + 29*IntegerPow<2>(z)*(210 - 868*IntegerPow<2>(z) + 1023*IntegerPow<4>(z)))))/32768.;
    if constexpr(I == 298) return (-15*sqrt(676039/pi)*y*(-IntegerPow<7>(x) + 7*IntegerPow<5>(x)*IntegerPow<2>(y) - 7*IntegerPow<3>(x)*IntegerPow<4>(y) + x*IntegerPow<6>(y))*z*(35 - 1260*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(126 - 372*IntegerPow<2>(z) + 341*IntegerPow<4>(z))))/4096.;
    if constexpr(I == 299) return (-3*sqrt(3380195/(2.*pi))*y*(-7*IntegerPow<6>(x) + 35*IntegerPow<4>(x)*IntegerPow<2>(y) - 21*IntegerPow<2>(x)*IntegerPow<4>(y) + IntegerPow<6>(y))*(-7 + 5*IntegerPow<2>(z)*(175 - 3150*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(210 - 465*IntegerPow<2>(z) + 341*IntegerPow<4>(z)))))/32768.;
    if constexpr(I == 300) return (sqrt(111546435/pi)*y*(3*IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 3*x*IntegerPow<4>(y))*z*(-21 + 875*IntegerPow<2>(z) - 9450*IntegerPow<4>(z) + 435*IntegerPow<6>(z)*(90 - 155*IntegerPow<2>(z) + 93*IntegerPow<4>(z))))/8192.;
    if constexpr(I == 301) return (3*sqrt(1616615/pi)*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(7 + 23*IntegerPow<2>(z)*(-42 + 875*IntegerPow<2>(z) - 6300*IntegerPow<4>(z) + 435*IntegerPow<6>(z)*(45 + 31*IntegerPow<2>(z)*(-2 + IntegerPow<2>(z))))))/32768.;
    if constexpr(I == 302) return (-3*sqrt(11305/(2.*pi))*x*y*z*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(1001 + 23*IntegerPow<2>(z)*(-2002 + 5*IntegerPow<2>(z)*(5005 - 25740*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(715 - 806*IntegerPow<2>(z) + 341*IntegerPow<4>(z))))))/4096.;
    if constexpr(I == 303) return -(sqrt(33915/pi)*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*(-143 + 21021*IntegerPow<2>(z) + 23*IntegerPow<4>(z)*(-21021 + 5*IntegerPow<2>(z)*(35035 + 3*IntegerPow<2>(z)*(-45045 + 87087*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(-91 + 33*IntegerPow<2>(z)))))))/32768.;
    if constexpr(I == 304) return (3*sqrt(11305/pi)*x*y*z*(-715 + 35035*IntegerPow<2>(z) + 23*IntegerPow<4>(z)*(-21021 + 5*IntegerPow<2>(z)*(25025 - 75075*IntegerPow<2>(z) + 118755*IntegerPow<4>(z) - 94395*IntegerPow<6>(z) + 29667*IntegerPow<8>(z)))))/8192.;
    if constexpr(I == 305) return (3*sqrt(595/pi)*y*(715 + 19*IntegerPow<2>(z)*(-5720 + 140140*IntegerPow<2>(z) + 23*IntegerPow<4>(z)*(-56056 + 5*IntegerPow<2>(z)*(50050 + 3*IntegerPow<2>(z)*(-40040 + 52780*IntegerPow<2>(z) - 35960*IntegerPow<4>(z) + 9889*IntegerPow<6>(z)))))))/65536.;
    if constexpr(I == 306) return (sqrt(35/pi)*z*(109395 + 19*IntegerPow<2>(z)*(-291720 + 4288284*IntegerPow<2>(z) + 23*IntegerPow<4>(z)*(-1225224 + 5*IntegerPow<2>(z)*(850850 + 9*IntegerPow<2>(z)*(-185640 + 29*IntegerPow<2>(z)*(7140 - 4216*IntegerPow<2>(z) + 1023*IntegerPow<4>(z))))))))/65536.;
    if constexpr(I == 307) return (3*sqrt(595/pi)*x*(715 + 19*IntegerPow<2>(z)*(-5720 + 140140*IntegerPow<2>(z) + 23*IntegerPow<4>(z)*(-56056 + 5*IntegerPow<2>(z)*(50050 + 3*IntegerPow<2>(z)*(-40040 + 52780*IntegerPow<2>(z) - 35960*IntegerPow<4>(z) + 9889*IntegerPow<6>(z)))))))/65536.;
    if constexpr(I == 308) return (3*sqrt(11305/pi)*(x - y)*(x + y)*z*(-715 + 35035*IntegerPow<2>(z) + 23*IntegerPow<4>(z)*(-21021 + 5*IntegerPow<2>(z)*(25025 - 75075*IntegerPow<2>(z) + 118755*IntegerPow<4>(z) - 94395*IntegerPow<6>(z) + 29667*IntegerPow<8>(z)))))/16384.;
    if constexpr(I == 309) return (sqrt(33915/pi)*(IntegerPow<3>(x) - 3*x*IntegerPow<2>(y))*(-143 + 21021*IntegerPow<2>(z) + 23*IntegerPow<4>(z)*(-21021 + 5*IntegerPow<2>(z)*(35035 + 3*IntegerPow<2>(z)*(-45045 + 87087*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(-91 + 33*IntegerPow<2>(z)))))))/32768.;
    if constexpr(I == 310) return (3*sqrt(11305/(2.*pi))*z*(8*IntegerPow<4>(y) + 8*IntegerPow<2>(y)*(-1 + IntegerPow<2>(z)) + IntegerPow<2>(-1 + IntegerPow<2>(z)))*(1001 + 23*IntegerPow<2>(z)*(-2002 + 5*IntegerPow<2>(z)*(5005 - 25740*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(715 - 806*IntegerPow<2>(z) + 341*IntegerPow<4>(z))))))/16384.;
    if constexpr(I == 311) return (3*sqrt(1616615/pi)*(IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 5*x*IntegerPow<4>(y))*(7 + 23*IntegerPow<2>(z)*(-42 + 875*IntegerPow<2>(z) - 6300*IntegerPow<4>(z) + 435*IntegerPow<6>(z)*(45 + 31*IntegerPow<2>(z)*(-2 + IntegerPow<2>(z))))))/32768.;
    if constexpr(I == 312) return -(sqrt(111546435/pi)*z*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(16*IntegerPow<4>(y) + 16*IntegerPow<2>(y)*(-1 + IntegerPow<2>(z)) + IntegerPow<2>(-1 + IntegerPow<2>(z)))*(-21 + 875*IntegerPow<2>(z) - 9450*IntegerPow<4>(z) + 435*IntegerPow<6>(z)*(90 - 155*IntegerPow<2>(z) + 93*IntegerPow<4>(z))))/16384.;
    if constexpr(I == 313) return (3*sqrt(3380195/(2.*pi))*(IntegerPow<7>(x) - 21*IntegerPow<5>(x)*IntegerPow<2>(y) + 35*IntegerPow<3>(x)*IntegerPow<4>(y) - 7*x*IntegerPow<6>(y))*(-7 + 5*IntegerPow<2>(z)*(175 - 3150*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(210 - 465*IntegerPow<2>(z) + 341*IntegerPow<4>(z)))))/32768.;
    if constexpr(I == 314) return (15*sqrt(676039/pi)*z*(128*IntegerPow<8>(y) + 256*IntegerPow<6>(y)*(-1 + IntegerPow<2>(z)) + 160*IntegerPow<4>(y)*IntegerPow<2>(-1 + IntegerPow<2>(z)) + 32*IntegerPow<2>(y)*IntegerPow<3>(-1 + IntegerPow<2>(z)) + IntegerPow<4>(-1 + IntegerPow<2>(z)))*(35 - 1260*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(126 - 372*IntegerPow<2>(z) + 341*IntegerPow<4>(z))))/32768.;
    if constexpr(I == 315) return (5*sqrt(52003/(2.*pi))*(IntegerPow<9>(x) - 36*IntegerPow<7>(x)*IntegerPow<2>(y) + 126*IntegerPow<5>(x)*IntegerPow<4>(y) - 84*IntegerPow<3>(x)*IntegerPow<6>(y) + 9*x*IntegerPow<8>(y))*(35 + 9*IntegerPow<2>(z)*(-420 + 29*IntegerPow<2>(z)*(210 - 868*IntegerPow<2>(z) + 1023*IntegerPow<4>(z)))))/32768.;
    if constexpr(I == 316) return (15*sqrt(156009/pi)*(IntegerPow<10>(x) - 45*IntegerPow<8>(x)*IntegerPow<2>(y) + 210*IntegerPow<6>(x)*IntegerPow<4>(y) - 210*IntegerPow<4>(x)*IntegerPow<6>(y) + 45*IntegerPow<2>(x)*IntegerPow<8>(y) - IntegerPow<10>(y))*z*(-35 + 29*IntegerPow<2>(z)*(35 - 217*IntegerPow<2>(z) + 341*IntegerPow<4>(z))))/16384.;
    if constexpr(I == 317) return (15*sqrt(156009/pi)*(IntegerPow<11>(x) - 55*IntegerPow<9>(x)*IntegerPow<2>(y) + 330*IntegerPow<7>(x)*IntegerPow<4>(y) - 462*IntegerPow<5>(x)*IntegerPow<6>(y) + 165*IntegerPow<3>(x)*IntegerPow<8>(y) - 11*x*IntegerPow<10>(y))*(-5 + 435*IntegerPow<2>(z) - 4495*IntegerPow<4>(z) + 9889*IntegerPow<6>(z)))/32768.;
    if constexpr(I == 318) return (15*sqrt(1508087/(2.*pi))*z*(15 - 310*IntegerPow<2>(z) + 1023*IntegerPow<4>(z))*(2048*IntegerPow<12>(y) + 6144*IntegerPow<10>(y)*(-1 + IntegerPow<2>(z)) + 6912*IntegerPow<8>(y)*IntegerPow<2>(-1 + IntegerPow<2>(z)) + 3584*IntegerPow<6>(y)*IntegerPow<3>(-1 + IntegerPow<2>(z)) + 840*IntegerPow<4>(y)*IntegerPow<4>(-1 + IntegerPow<2>(z)) + 72*IntegerPow<2>(y)*IntegerPow<5>(-1 + IntegerPow<2>(z)) + IntegerPow<6>(-1 + IntegerPow<2>(z))))/16384.;
    if constexpr(I == 319) return (15*sqrt(4524261/pi)*(IntegerPow<13>(x) - 78*IntegerPow<11>(x)*IntegerPow<2>(y) + 715*IntegerPow<9>(x)*IntegerPow<4>(y) - 1716*IntegerPow<7>(x)*IntegerPow<6>(y) + 1287*IntegerPow<5>(x)*IntegerPow<8>(y) - 286*IntegerPow<3>(x)*IntegerPow<10>(y) + 13*x*IntegerPow<12>(y))*(1 - 62*IntegerPow<2>(z) + 341*IntegerPow<4>(z)))/32768.;
    if constexpr(I == 320) return (15*sqrt(140252091/pi)*(IntegerPow<14>(x) - 91*IntegerPow<12>(x)*IntegerPow<2>(y) + 1001*IntegerPow<10>(x)*IntegerPow<4>(y) - 3003*IntegerPow<8>(x)*IntegerPow<6>(y) + 3003*IntegerPow<6>(x)*IntegerPow<8>(y) - 1001*IntegerPow<4>(x)*IntegerPow<10>(y) + 91*IntegerPow<2>(x)*IntegerPow<12>(y) - IntegerPow<14>(y))*z*(-1 + 11*IntegerPow<2>(z)))/16384.;
    if constexpr(I == 321) return (15*sqrt(46750697/(2.*pi))*x*(IntegerPow<2>(x) - 3*IntegerPow<2>(y))*(IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + 5*IntegerPow<4>(y))*(IntegerPow<8>(x) - 92*IntegerPow<6>(x)*IntegerPow<2>(y) + 134*IntegerPow<4>(x)*IntegerPow<4>(y) - 28*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*(-1 + 33*IntegerPow<2>(z)))/65536.;
    if constexpr(I == 322) return (15*sqrt(1542773001/pi)*(IntegerPow<16>(x) - 120*IntegerPow<14>(x)*IntegerPow<2>(y) + 1820*IntegerPow<12>(x)*IntegerPow<4>(y) - 8008*IntegerPow<10>(x)*IntegerPow<6>(y) + 12870*IntegerPow<8>(x)*IntegerPow<8>(y) - 8008*IntegerPow<6>(x)*IntegerPow<10>(y) + 1820*IntegerPow<4>(x)*IntegerPow<12>(y) - 120*IntegerPow<2>(x)*IntegerPow<14>(y) + IntegerPow<16>(y))*z)/65536.;
    if constexpr(I == 323) return (15*sqrt(90751353/(2.*pi))*(IntegerPow<17>(x) - 136*IntegerPow<15>(x)*IntegerPow<2>(y) + 2380*IntegerPow<13>(x)*IntegerPow<4>(y) - 12376*IntegerPow<11>(x)*IntegerPow<6>(y) + 24310*IntegerPow<9>(x)*IntegerPow<8>(y) - 19448*IntegerPow<7>(x)*IntegerPow<10>(y) + 6188*IntegerPow<5>(x)*IntegerPow<12>(y) - 680*IntegerPow<3>(x)*IntegerPow<14>(y) + 17*x*IntegerPow<16>(y)))/65536.;
    if constexpr(I == 324) return (5*sqrt(3357800061/(2.*pi))*x*y*(IntegerPow<2>(x) - 3*IntegerPow<2>(y))*(3*IntegerPow<2>(x) - IntegerPow<2>(y))*(IntegerPow<6>(x) - 33*IntegerPow<4>(x)*IntegerPow<2>(y) + 27*IntegerPow<2>(x)*IntegerPow<4>(y) - 3*IntegerPow<6>(y))*(3*IntegerPow<6>(x) - 27*IntegerPow<4>(x)*IntegerPow<2>(y) + 33*IntegerPow<2>(x)*IntegerPow<4>(y) - IntegerPow<6>(y)))/65536.;
    if constexpr(I == 325) return (15*sqrt(3357800061/(2.*pi))*y*(17*IntegerPow<16>(x) - 680*IntegerPow<14>(x)*IntegerPow<2>(y) + 6188*IntegerPow<12>(x)*IntegerPow<4>(y) - 19448*IntegerPow<10>(x)*IntegerPow<6>(y) + 24310*IntegerPow<8>(x)*IntegerPow<8>(y) - 12376*IntegerPow<6>(x)*IntegerPow<10>(y) + 2380*IntegerPow<4>(x)*IntegerPow<12>(y) - 136*IntegerPow<2>(x)*IntegerPow<14>(y) + IntegerPow<16>(y))*z)/65536.;
    if constexpr(I == 326) return (-3*sqrt(2398428615/pi)*y*(-IntegerPow<15>(x) + 35*IntegerPow<13>(x)*IntegerPow<2>(y) - 273*IntegerPow<11>(x)*IntegerPow<4>(y) + 715*IntegerPow<9>(x)*IntegerPow<6>(y) - 715*IntegerPow<7>(x)*IntegerPow<8>(y) + 273*IntegerPow<5>(x)*IntegerPow<10>(y) - 35*IntegerPow<3>(x)*IntegerPow<12>(y) + x*IntegerPow<14>(y))*(-1 + 35*IntegerPow<2>(z)))/8192.;
    if constexpr(I == 327) return (-3*sqrt(13591095485/(2.*pi))*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(IntegerPow<8>(x) - 28*IntegerPow<6>(x)*IntegerPow<2>(y) + 134*IntegerPow<4>(x)*IntegerPow<4>(y) - 92*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*z*(-3 + 35*IntegerPow<2>(z)))/65536.;
    if constexpr(I == 328) return (3*sqrt(3706662405/(2.*pi))*x*y*(IntegerPow<6>(x) - 21*IntegerPow<4>(x)*IntegerPow<2>(y) + 35*IntegerPow<2>(x)*IntegerPow<4>(y) - 7*IntegerPow<6>(y))*(7*IntegerPow<6>(x) - 35*IntegerPow<4>(x)*IntegerPow<2>(y) + 21*IntegerPow<2>(x)*IntegerPow<4>(y) - IntegerPow<6>(y))*(1 - 66*IntegerPow<2>(z) + 385*IntegerPow<4>(z)))/65536.;
    if constexpr(I == 329) return (15*sqrt(741332481/pi)*y*(13*IntegerPow<12>(x) - 286*IntegerPow<10>(x)*IntegerPow<2>(y) + 1287*IntegerPow<8>(x)*IntegerPow<4>(y) - 1716*IntegerPow<6>(x)*IntegerPow<6>(y) + 715*IntegerPow<4>(x)*IntegerPow<8>(y) - 78*IntegerPow<2>(x)*IntegerPow<10>(y) + IntegerPow<12>(y))*(z - 22*IntegerPow<3>(z) + 77*IntegerPow<5>(z)))/32768.;
    if constexpr(I == 330) return (15*sqrt(7971317/(2.*pi))*x*(x - y)*y*(x + y)*(IntegerPow<2>(x) - 3*IntegerPow<2>(y))*(3*IntegerPow<2>(x) - IntegerPow<2>(y))*(IntegerPow<2>(x) - 4*x*y + IntegerPow<2>(y))*(IntegerPow<2>(x) + 4*x*y + IntegerPow<2>(y))*(-1 + 31*IntegerPow<2>(z)*(3 - 33*IntegerPow<2>(z) + 77*IntegerPow<4>(z))))/8192.;
    if constexpr(I == 331) return (-3*sqrt(836988285/pi)*y*(-11*IntegerPow<10>(x) + 165*IntegerPow<8>(x)*IntegerPow<2>(y) - 462*IntegerPow<6>(x)*IntegerPow<4>(y) + 330*IntegerPow<4>(x)*IntegerPow<6>(y) - 55*IntegerPow<2>(x)*IntegerPow<8>(y) + IntegerPow<10>(y))*z*(-5 + 31*IntegerPow<2>(z)*(5 - 33*IntegerPow<2>(z) + 55*IntegerPow<4>(z))))/32768.;
    if constexpr(I == 332) return (3*sqrt(28861665/(2.*pi))*x*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + 5*IntegerPow<4>(y))*(5 - 580*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(10 - 44*IntegerPow<2>(z) + 55*IntegerPow<4>(z))))/32768.;
    if constexpr(I == 333) return (sqrt(4123095/(2.*pi))*y*(9*IntegerPow<8>(x) - 84*IntegerPow<6>(x)*IntegerPow<2>(y) + 126*IntegerPow<4>(x)*IntegerPow<4>(y) - 36*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*z*(315 + 29*IntegerPow<2>(z)*(-420 + 31*IntegerPow<2>(z)*(126 - 396*IntegerPow<2>(z) + 385*IntegerPow<4>(z)))))/32768.;
    if constexpr(I == 334) return (-15*sqrt(274873/pi)*y*(-IntegerPow<7>(x) + 7*IntegerPow<5>(x)*IntegerPow<2>(y) - 7*IntegerPow<3>(x)*IntegerPow<4>(y) + x*IntegerPow<6>(y))*(-7 + 3*IntegerPow<2>(z)*(315 - 6090*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(42 - 99*IntegerPow<2>(z) + 77*IntegerPow<4>(z)))))/8192.;
    if constexpr(I == 335) return (-15*sqrt(39306839/(2.*pi))*y*(-7*IntegerPow<6>(x) + 35*IntegerPow<4>(x)*IntegerPow<2>(y) - 21*IntegerPow<2>(x)*IntegerPow<4>(y) + IntegerPow<6>(y))*z*(-7 + 315*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(-42 + 31*IntegerPow<2>(z)*(6 - 11*IntegerPow<2>(z) + 7*IntegerPow<4>(z)))))/32768.;
    if constexpr(I == 336) return (sqrt(117920517/(2.*pi))*y*(3*IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 3*x*IntegerPow<4>(y))*(7 + 15*IntegerPow<2>(z)*(-70 + 1575*IntegerPow<2>(z) + 29*IntegerPow<4>(z)*(-420 + 31*IntegerPow<2>(z)*(45 - 66*IntegerPow<2>(z) + 35*IntegerPow<4>(z))))))/32768.;
    if constexpr(I == 337) return (3*sqrt(3023603/pi)*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*z*(91 + 5*IntegerPow<2>(z)*(-910 + 12285*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(-780 + 31*IntegerPow<2>(z)*(65 - 78*IntegerPow<2>(z) + 35*IntegerPow<4>(z))))))/32768.;
    if constexpr(I == 338) return (-3*sqrt(920227/(2.*pi))*x*y*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(-13 + 23*IntegerPow<2>(z)*(91 + 5*IntegerPow<2>(z)*(-455 + 3*IntegerPow<2>(z)*(1365 - 5655*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(13 - 13*IntegerPow<2>(z) + 5*IntegerPow<4>(z)))))))/8192.;
    if constexpr(I == 339) return -(sqrt(1254855/pi)*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*z*(-429 + 23023*IntegerPow<2>(z) + 345*IntegerPow<4>(z)*(-1001 + 6435*IntegerPow<2>(z) + 29*IntegerPow<4>(z)*(-715 + 31*IntegerPow<2>(z)*(39 + 11*IntegerPow<2>(z)*(-3 + IntegerPow<2>(z)))))))/32768.;
    if constexpr(I == 340) return (3*sqrt(59755/pi)*x*y*(143 - 24024*IntegerPow<2>(z) + 23*IntegerPow<4>(z)*(28028 + 5*IntegerPow<2>(z)*(-56056 + 270270*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(-8008 + 31*IntegerPow<2>(z)*(364 - 264*IntegerPow<2>(z) + 77*IntegerPow<4>(z)))))))/65536.;
    if constexpr(I == 341) return (3*sqrt(703/pi)*y*z*(12155 - 680680*IntegerPow<2>(z) + 23*IntegerPow<4>(z)*(476476 + 5*IntegerPow<2>(z)*(-680680 + 2552550*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(-61880 + 31*IntegerPow<2>(z)*(2380 - 1496*IntegerPow<2>(z) + 385*IntegerPow<4>(z)))))))/65536.;
    if constexpr(I == 342) return (sqrt(37/pi)*(-12155 + 57*IntegerPow<2>(z)*(36465 - 1021020*IntegerPow<2>(z) + 23*IntegerPow<4>(z)*(476476 + 5*IntegerPow<2>(z)*(-510510 + 1531530*IntegerPow<2>(z) + 29*IntegerPow<4>(z)*(-92820 + 31*IntegerPow<2>(z)*(3060 - 1683*IntegerPow<2>(z) + 385*IntegerPow<4>(z))))))))/131072.;
    if constexpr(I == 343) return (3*sqrt(703/pi)*x*z*(12155 - 680680*IntegerPow<2>(z) + 23*IntegerPow<4>(z)*(476476 + 5*IntegerPow<2>(z)*(-680680 + 2552550*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(-61880 + 31*IntegerPow<2>(z)*(2380 - 1496*IntegerPow<2>(z) + 385*IntegerPow<4>(z)))))))/65536.;
    if constexpr(I == 344) return (-3*sqrt(59755/pi)*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(143 - 24024*IntegerPow<2>(z) + 23*IntegerPow<4>(z)*(28028 + 5*IntegerPow<2>(z)*(-56056 + 270270*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(-8008 + 31*IntegerPow<2>(z)*(364 - 264*IntegerPow<2>(z) + 77*IntegerPow<4>(z)))))))/131072.;
    if constexpr(I == 345) return (sqrt(1254855/pi)*(IntegerPow<3>(x) - 3*x*IntegerPow<2>(y))*z*(-429 + 23023*IntegerPow<2>(z) + 345*IntegerPow<4>(z)*(-1001 + 6435*IntegerPow<2>(z) + 29*IntegerPow<4>(z)*(-715 + 31*IntegerPow<2>(z)*(39 + 11*IntegerPow<2>(z)*(-3 + IntegerPow<2>(z)))))))/32768.;
    if constexpr(I == 346) return (3*sqrt(920227/(2.*pi))*(IntegerPow<4>(x) - 6*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(-13 + 23*IntegerPow<2>(z)*(91 + 5*IntegerPow<2>(z)*(-455 + 3*IntegerPow<2>(z)*(1365 - 5655*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(13 - 13*IntegerPow<2>(z) + 5*IntegerPow<4>(z)))))))/32768.;
    if constexpr(I == 347) return (3*sqrt(3023603/pi)*(IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 5*x*IntegerPow<4>(y))*z*(91 + 5*IntegerPow<2>(z)*(-910 + 12285*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(-780 + 31*IntegerPow<2>(z)*(65 - 78*IntegerPow<2>(z) + 35*IntegerPow<4>(z))))))/32768.;
    if constexpr(I == 348) return -(sqrt(117920517/(2.*pi))*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(16*IntegerPow<4>(y) + 16*IntegerPow<2>(y)*(-1 + IntegerPow<2>(z)) + IntegerPow<2>(-1 + IntegerPow<2>(z)))*(7 + 15*IntegerPow<2>(z)*(-70 + 1575*IntegerPow<2>(z) + 29*IntegerPow<4>(z)*(-420 + 31*IntegerPow<2>(z)*(45 - 66*IntegerPow<2>(z) + 35*IntegerPow<4>(z))))))/65536.;
    if constexpr(I == 349) return (15*sqrt(39306839/(2.*pi))*(IntegerPow<7>(x) - 21*IntegerPow<5>(x)*IntegerPow<2>(y) + 35*IntegerPow<3>(x)*IntegerPow<4>(y) - 7*x*IntegerPow<6>(y))*z*(-7 + 315*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(-42 + 31*IntegerPow<2>(z)*(6 - 11*IntegerPow<2>(z) + 7*IntegerPow<4>(z)))))/32768.;
    if constexpr(I == 350) return (15*sqrt(274873/pi)*(128*IntegerPow<8>(y) + 256*IntegerPow<6>(y)*(-1 + IntegerPow<2>(z)) + 160*IntegerPow<4>(y)*IntegerPow<2>(-1 + IntegerPow<2>(z)) + 32*IntegerPow<2>(y)*IntegerPow<3>(-1 + IntegerPow<2>(z)) + IntegerPow<4>(-1 + IntegerPow<2>(z)))*(-7 + 3*IntegerPow<2>(z)*(315 - 6090*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(42 - 99*IntegerPow<2>(z) + 77*IntegerPow<4>(z)))))/65536.;
    if constexpr(I == 351) return (sqrt(4123095/(2.*pi))*(IntegerPow<9>(x) - 36*IntegerPow<7>(x)*IntegerPow<2>(y) + 126*IntegerPow<5>(x)*IntegerPow<4>(y) - 84*IntegerPow<3>(x)*IntegerPow<6>(y) + 9*x*IntegerPow<8>(y))*z*(315 + 29*IntegerPow<2>(z)*(-420 + 31*IntegerPow<2>(z)*(126 - 396*IntegerPow<2>(z) + 385*IntegerPow<4>(z)))))/32768.;
    if constexpr(I == 352) return (-3*sqrt(28861665/(2.*pi))*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(256*IntegerPow<8>(y) + 512*IntegerPow<6>(y)*(-1 + IntegerPow<2>(z)) + 304*IntegerPow<4>(y)*IntegerPow<2>(-1 + IntegerPow<2>(z)) + 48*IntegerPow<2>(y)*IntegerPow<3>(-1 + IntegerPow<2>(z)) + IntegerPow<4>(-1 + IntegerPow<2>(z)))*(5 - 580*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(10 - 44*IntegerPow<2>(z) + 55*IntegerPow<4>(z))))/65536.;
    if constexpr(I == 353) return (3*sqrt(836988285/pi)*(IntegerPow<11>(x) - 55*IntegerPow<9>(x)*IntegerPow<2>(y) + 330*IntegerPow<7>(x)*IntegerPow<4>(y) - 462*IntegerPow<5>(x)*IntegerPow<6>(y) + 165*IntegerPow<3>(x)*IntegerPow<8>(y) - 11*x*IntegerPow<10>(y))*z*(-5 + 31*IntegerPow<2>(z)*(5 - 33*IntegerPow<2>(z) + 55*IntegerPow<4>(z))))/32768.;
    if constexpr(I == 354) return (15*sqrt(7971317/(2.*pi))*(IntegerPow<12>(x) - 66*IntegerPow<10>(x)*IntegerPow<2>(y) + 495*IntegerPow<8>(x)*IntegerPow<4>(y) - 924*IntegerPow<6>(x)*IntegerPow<6>(y) + 495*IntegerPow<4>(x)*IntegerPow<8>(y) - 66*IntegerPow<2>(x)*IntegerPow<10>(y) + IntegerPow<12>(y))*(-1 + 31*IntegerPow<2>(z)*(3 - 33*IntegerPow<2>(z) + 77*IntegerPow<4>(z))))/32768.;
    if constexpr(I == 355) return (15*sqrt(741332481/pi)*(IntegerPow<13>(x) - 78*IntegerPow<11>(x)*IntegerPow<2>(y) + 715*IntegerPow<9>(x)*IntegerPow<4>(y) - 1716*IntegerPow<7>(x)*IntegerPow<6>(y) + 1287*IntegerPow<5>(x)*IntegerPow<8>(y) - 286*IntegerPow<3>(x)*IntegerPow<10>(y) + 13*x*IntegerPow<12>(y))*(z - 22*IntegerPow<3>(z) + 77*IntegerPow<5>(z)))/32768.;
    if constexpr(I == 356) return (3*sqrt(3706662405/(2.*pi))*(IntegerPow<14>(x) - 91*IntegerPow<12>(x)*IntegerPow<2>(y) + 1001*IntegerPow<10>(x)*IntegerPow<4>(y) - 3003*IntegerPow<8>(x)*IntegerPow<6>(y) + 3003*IntegerPow<6>(x)*IntegerPow<8>(y) - 1001*IntegerPow<4>(x)*IntegerPow<10>(y) + 91*IntegerPow<2>(x)*IntegerPow<12>(y) - IntegerPow<14>(y))*(1 - 66*IntegerPow<2>(z) + 385*IntegerPow<4>(z)))/131072.;
    if constexpr(I == 357) return (3*sqrt(13591095485/(2.*pi))*x*(IntegerPow<2>(x) - 3*IntegerPow<2>(y))*(IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + 5*IntegerPow<4>(y))*(IntegerPow<8>(x) - 92*IntegerPow<6>(x)*IntegerPow<2>(y) + 134*IntegerPow<4>(x)*IntegerPow<4>(y) - 28*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*z*(-3 + 35*IntegerPow<2>(z)))/65536.;
    if constexpr(I == 358) return (3*sqrt(2398428615/pi)*(IntegerPow<16>(x) - 120*IntegerPow<14>(x)*IntegerPow<2>(y) + 1820*IntegerPow<12>(x)*IntegerPow<4>(y) - 8008*IntegerPow<10>(x)*IntegerPow<6>(y) + 12870*IntegerPow<8>(x)*IntegerPow<8>(y) - 8008*IntegerPow<6>(x)*IntegerPow<10>(y) + 1820*IntegerPow<4>(x)*IntegerPow<12>(y) - 120*IntegerPow<2>(x)*IntegerPow<14>(y) + IntegerPow<16>(y))*(-1 + 35*IntegerPow<2>(z)))/131072.;
    if constexpr(I == 359) return (15*sqrt(3357800061/(2.*pi))*(IntegerPow<17>(x) - 136*IntegerPow<15>(x)*IntegerPow<2>(y) + 2380*IntegerPow<13>(x)*IntegerPow<4>(y) - 12376*IntegerPow<11>(x)*IntegerPow<6>(y) + 24310*IntegerPow<9>(x)*IntegerPow<8>(y) - 19448*IntegerPow<7>(x)*IntegerPow<10>(y) + 6188*IntegerPow<5>(x)*IntegerPow<12>(y) - 680*IntegerPow<3>(x)*IntegerPow<14>(y) + 17*x*IntegerPow<16>(y))*z)/65536.;
    if constexpr(I == 360) return (5*sqrt(3357800061/(2.*pi))*(IntegerPow<18>(x) - 153*IntegerPow<16>(x)*IntegerPow<2>(y) + 3060*IntegerPow<14>(x)*IntegerPow<4>(y) - 18564*IntegerPow<12>(x)*IntegerPow<6>(y) + 43758*IntegerPow<10>(x)*IntegerPow<8>(y) - 43758*IntegerPow<8>(x)*IntegerPow<10>(y) + 18564*IntegerPow<6>(x)*IntegerPow<12>(y) - 3060*IntegerPow<4>(x)*IntegerPow<14>(y) + 153*IntegerPow<2>(x)*IntegerPow<16>(y) - IntegerPow<18>(y)))/131072.;
    if constexpr(I == 361) return (-15*sqrt(765814049/pi)*y*(-19*IntegerPow<18>(x) + 969*IntegerPow<16>(x)*IntegerPow<2>(y) - 11628*IntegerPow<14>(x)*IntegerPow<4>(y) + 50388*IntegerPow<12>(x)*IntegerPow<6>(y) - 92378*IntegerPow<10>(x)*IntegerPow<8>(y) + 75582*IntegerPow<8>(x)*IntegerPow<10>(y) - 27132*IntegerPow<6>(x)*IntegerPow<12>(y) + 3876*IntegerPow<4>(x)*IntegerPow<14>(y) - 171*IntegerPow<2>(x)*IntegerPow<16>(y) + IntegerPow<18>(y)))/262144.;
    if constexpr(I == 362) return (15*sqrt(14550466931/(2.*pi))*x*y*(IntegerPow<2>(x) - 3*IntegerPow<2>(y))*(3*IntegerPow<2>(x) - IntegerPow<2>(y))*(IntegerPow<6>(x) - 33*IntegerPow<4>(x)*IntegerPow<2>(y) + 27*IntegerPow<2>(x)*IntegerPow<4>(y) - 3*IntegerPow<6>(y))*(3*IntegerPow<6>(x) - 27*IntegerPow<4>(x)*IntegerPow<2>(y) + 33*IntegerPow<2>(x)*IntegerPow<4>(y) - IntegerPow<6>(y))*z)/65536.;
    if constexpr(I == 363) return (15*sqrt(393255863/pi)*y*(17*IntegerPow<16>(x) - 680*IntegerPow<14>(x)*IntegerPow<2>(y) + 6188*IntegerPow<12>(x)*IntegerPow<4>(y) - 19448*IntegerPow<10>(x)*IntegerPow<6>(y) + 24310*IntegerPow<8>(x)*IntegerPow<8>(y) - 12376*IntegerPow<6>(x)*IntegerPow<10>(y) + 2380*IntegerPow<4>(x)*IntegerPow<12>(y) - 136*IntegerPow<2>(x)*IntegerPow<14>(y) + IntegerPow<16>(y))*(-1 + 37*IntegerPow<2>(z)))/262144.;
    if constexpr(I == 364) return (-15*sqrt(1179767589/pi)*y*(-IntegerPow<15>(x) + 35*IntegerPow<13>(x)*IntegerPow<2>(y) - 273*IntegerPow<11>(x)*IntegerPow<4>(y) + 715*IntegerPow<9>(x)*IntegerPow<6>(y) - 715*IntegerPow<7>(x)*IntegerPow<8>(y) + 273*IntegerPow<5>(x)*IntegerPow<10>(y) - 35*IntegerPow<3>(x)*IntegerPow<12>(y) + x*IntegerPow<14>(y))*z*(-3 + 37*IntegerPow<2>(z)))/8192.;
    if constexpr(I == 365) return (-3*sqrt(842691135/pi)*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(IntegerPow<8>(x) - 28*IntegerPow<6>(x)*IntegerPow<2>(y) + 134*IntegerPow<4>(x)*IntegerPow<4>(y) - 92*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*(3 + 35*IntegerPow<2>(z)*(-6 + 37*IntegerPow<2>(z))))/262144.;
    if constexpr(I == 366) return (15*sqrt(2865149859/(2.*pi))*x*y*(IntegerPow<6>(x) - 21*IntegerPow<4>(x)*IntegerPow<2>(y) + 35*IntegerPow<2>(x)*IntegerPow<4>(y) - 7*IntegerPow<6>(y))*(7*IntegerPow<6>(x) - 35*IntegerPow<4>(x)*IntegerPow<2>(y) + 21*IntegerPow<2>(x)*IntegerPow<4>(y) - IntegerPow<6>(y))*z*(3 - 70*IntegerPow<2>(z) + 259*IntegerPow<4>(z)))/65536.;
    if constexpr(I == 367) return (15*sqrt(260468169/pi)*y*(13*IntegerPow<12>(x) - 286*IntegerPow<10>(x)*IntegerPow<2>(y) + 1287*IntegerPow<8>(x)*IntegerPow<4>(y) - 1716*IntegerPow<6>(x)*IntegerPow<6>(y) + 715*IntegerPow<4>(x)*IntegerPow<8>(y) - 78*IntegerPow<2>(x)*IntegerPow<10>(y) + IntegerPow<12>(y))*(-1 + 99*IntegerPow<2>(z) - 1155*IntegerPow<4>(z) + 2849*IntegerPow<6>(z)))/262144.;
    if constexpr(I == 368) return (15*sqrt(1823277183/(2.*pi))*x*(x - y)*y*(x + y)*(IntegerPow<2>(x) - 3*IntegerPow<2>(y))*(3*IntegerPow<2>(x) - IntegerPow<2>(y))*(IntegerPow<2>(x) - 4*x*y + IntegerPow<2>(y))*(IntegerPow<2>(x) + 4*x*y + IntegerPow<2>(y))*z*(-1 + 33*IntegerPow<2>(z) - 231*IntegerPow<4>(z) + 407*IntegerPow<6>(z)))/8192.;
    if constexpr(I == 369) return (-15*sqrt(58815393/pi)*y*(-11*IntegerPow<10>(x) + 165*IntegerPow<8>(x)*IntegerPow<2>(y) - 462*IntegerPow<6>(x)*IntegerPow<4>(y) + 330*IntegerPow<4>(x)*IntegerPow<6>(y) - 55*IntegerPow<2>(x)*IntegerPow<8>(y) + IntegerPow<10>(y))*(1 + 31*IntegerPow<2>(z)*(-4 + 66*IntegerPow<2>(z) - 308*IntegerPow<4>(z) + 407*IntegerPow<6>(z))))/131072.;
    if constexpr(I == 370) return (3*sqrt(98025655/(2.*pi))*x*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + 5*IntegerPow<4>(y))*z*(45 + 31*IntegerPow<2>(z)*(-60 + 594*IntegerPow<2>(z) - 1980*IntegerPow<4>(z) + 2035*IntegerPow<6>(z))))/32768.;
    if constexpr(I == 371) return (15*sqrt(676039/pi)*y*(9*IntegerPow<8>(x) - 84*IntegerPow<6>(x)*IntegerPow<2>(y) + 126*IntegerPow<4>(x)*IntegerPow<4>(y) - 36*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*(-9 + 29*IntegerPow<2>(z)*(45 + 31*IntegerPow<2>(z)*(-30 + 198*IntegerPow<2>(z) - 495*IntegerPow<4>(z) + 407*IntegerPow<6>(z)))))/131072.;
    if constexpr(I == 372) return (-15*sqrt(1062347/pi)*y*(-IntegerPow<7>(x) + 7*IntegerPow<5>(x)*IntegerPow<2>(y) - 7*IntegerPow<3>(x)*IntegerPow<4>(y) + x*IntegerPow<6>(y))*z*(-63 + 3045*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(-42 + 198*IntegerPow<2>(z) - 385*IntegerPow<4>(z) + 259*IntegerPow<6>(z))))/8192.;
    if constexpr(I == 373) return (-15*sqrt(1062347/pi)*y*(-7*IntegerPow<6>(x) + 35*IntegerPow<4>(x)*IntegerPow<2>(y) - 21*IntegerPow<2>(x)*IntegerPow<4>(y) + IntegerPow<6>(y))*(7 + 3*IntegerPow<2>(z)*(-378 + 9135*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(-84 + 297*IntegerPow<2>(z) - 462*IntegerPow<4>(z) + 259*IntegerPow<6>(z)))))/131072.;
    if constexpr(I == 374) return (15*sqrt(1062347/(2.*pi))*y*(3*IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 3*x*IntegerPow<4>(y))*z*(91 - 4914*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(819 + 31*IntegerPow<2>(z)*(-156 + 429*IntegerPow<2>(z) - 546*IntegerPow<4>(z) + 259*IntegerPow<6>(z)))))/32768.;
    if constexpr(I == 375) return (3*sqrt(7436429/pi)*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(-13 + 5*IntegerPow<2>(z)*(455 + 3*IntegerPow<2>(z)*(-4095 + 39585*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(-195 + 429*IntegerPow<2>(z) - 455*IntegerPow<4>(z) + 185*IntegerPow<6>(z))))))/131072.;
    if constexpr(I == 376) return (-3*sqrt(37182145/(2.*pi))*x*y*z*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(-39 + 5*IntegerPow<2>(z)*(455 - 7371*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(585 + 31*IntegerPow<2>(z)*(-65 + 117*IntegerPow<2>(z) - 105*IntegerPow<4>(z) + 37*IntegerPow<6>(z))))))/8192.;
    if constexpr(I == 377) return (-3*sqrt(1616615/(2.*pi))*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*(39 + 23*IntegerPow<2>(z)*(-312 + 9100*IntegerPow<2>(z) + 15*IntegerPow<4>(z)*(-6552 + 33930*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(-104 + 156*IntegerPow<2>(z) - 120*IntegerPow<4>(z) + 37*IntegerPow<6>(z))))))/131072.;
    if constexpr(I == 378) return (3*sqrt(8645/pi)*x*y*z*(7293 + 23*IntegerPow<2>(z)*(-19448 + 5*IntegerPow<2>(z)*(68068 + 3*IntegerPow<2>(z)*(-175032 + 704990*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(-1768 + 2244*IntegerPow<2>(z) - 1496*IntegerPow<4>(z) + 407*IntegerPow<6>(z)))))))/65536.;
    if constexpr(I == 379) return (sqrt(3705/(2.*pi))*y*(-2431 + 459459*IntegerPow<2>(z) + 69*IntegerPow<4>(z)*(-204204 + 5*IntegerPow<2>(z)*(476476 - 2756754*IntegerPow<2>(z) + 8882874*IntegerPow<4>(z) + 899*IntegerPow<6>(z)*(-18564 + 20196*IntegerPow<2>(z) - 11781*IntegerPow<4>(z) + 2849*IntegerPow<6>(z))))))/131072.;
    if constexpr(I == 380) return (sqrt(39/pi)*z*(-230945 + 14549535*IntegerPow<2>(z) + 69*IntegerPow<4>(z)*(-3879876 + 5*IntegerPow<2>(z)*(6466460 - 29099070*IntegerPow<2>(z) + 76715730*IntegerPow<4>(z) + 899*IntegerPow<6>(z)*(-135660 + 127908*IntegerPow<2>(z) - 65835*IntegerPow<4>(z) + 14245*IntegerPow<6>(z))))))/131072.;
    if constexpr(I == 381) return (sqrt(3705/(2.*pi))*x*(-2431 + 459459*IntegerPow<2>(z) + 69*IntegerPow<4>(z)*(-204204 + 5*IntegerPow<2>(z)*(476476 - 2756754*IntegerPow<2>(z) + 8882874*IntegerPow<4>(z) + 899*IntegerPow<6>(z)*(-18564 + 20196*IntegerPow<2>(z) - 11781*IntegerPow<4>(z) + 2849*IntegerPow<6>(z))))))/131072.;
    if constexpr(I == 382) return (-3*sqrt(8645/pi)*z*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(7293 + 23*IntegerPow<2>(z)*(-19448 + 5*IntegerPow<2>(z)*(68068 + 3*IntegerPow<2>(z)*(-175032 + 704990*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(-1768 + 2244*IntegerPow<2>(z) - 1496*IntegerPow<4>(z) + 407*IntegerPow<6>(z)))))))/131072.;
    if constexpr(I == 383) return (3*sqrt(1616615/(2.*pi))*(IntegerPow<3>(x) - 3*x*IntegerPow<2>(y))*(39 + 23*IntegerPow<2>(z)*(-312 + 9100*IntegerPow<2>(z) + 15*IntegerPow<4>(z)*(-6552 + 33930*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(-104 + 156*IntegerPow<2>(z) - 120*IntegerPow<4>(z) + 37*IntegerPow<6>(z))))))/131072.;
    if constexpr(I == 384) return (3*sqrt(37182145/(2.*pi))*(IntegerPow<4>(x) - 6*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*z*(-39 + 5*IntegerPow<2>(z)*(455 - 7371*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(585 + 31*IntegerPow<2>(z)*(-65 + 117*IntegerPow<2>(z) - 105*IntegerPow<4>(z) + 37*IntegerPow<6>(z))))))/32768.;
    if constexpr(I == 385) return (3*sqrt(7436429/pi)*(IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 5*x*IntegerPow<4>(y))*(-13 + 5*IntegerPow<2>(z)*(455 + 3*IntegerPow<2>(z)*(-4095 + 39585*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(-195 + 429*IntegerPow<2>(z) - 455*IntegerPow<4>(z) + 185*IntegerPow<6>(z))))))/131072.;
    if constexpr(I == 386) return (-15*sqrt(1062347/(2.*pi))*z*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(16*IntegerPow<4>(y) + 16*IntegerPow<2>(y)*(-1 + IntegerPow<2>(z)) + IntegerPow<2>(-1 + IntegerPow<2>(z)))*(91 - 4914*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(819 + 31*IntegerPow<2>(z)*(-156 + 429*IntegerPow<2>(z) - 546*IntegerPow<4>(z) + 259*IntegerPow<6>(z)))))/65536.;
    if constexpr(I == 387) return (15*sqrt(1062347/pi)*(IntegerPow<7>(x) - 21*IntegerPow<5>(x)*IntegerPow<2>(y) + 35*IntegerPow<3>(x)*IntegerPow<4>(y) - 7*x*IntegerPow<6>(y))*(7 + 3*IntegerPow<2>(z)*(-378 + 9135*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(-84 + 297*IntegerPow<2>(z) - 462*IntegerPow<4>(z) + 259*IntegerPow<6>(z)))))/131072.;
    if constexpr(I == 388) return (15*sqrt(1062347/pi)*z*(128*IntegerPow<8>(y) + 256*IntegerPow<6>(y)*(-1 + IntegerPow<2>(z)) + 160*IntegerPow<4>(y)*IntegerPow<2>(-1 + IntegerPow<2>(z)) + 32*IntegerPow<2>(y)*IntegerPow<3>(-1 + IntegerPow<2>(z)) + IntegerPow<4>(-1 + IntegerPow<2>(z)))*(-63 + 3045*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(-42 + 198*IntegerPow<2>(z) - 385*IntegerPow<4>(z) + 259*IntegerPow<6>(z))))/65536.;
    if constexpr(I == 389) return (15*sqrt(676039/pi)*(IntegerPow<9>(x) - 36*IntegerPow<7>(x)*IntegerPow<2>(y) + 126*IntegerPow<5>(x)*IntegerPow<4>(y) - 84*IntegerPow<3>(x)*IntegerPow<6>(y) + 9*x*IntegerPow<8>(y))*(-9 + 29*IntegerPow<2>(z)*(45 + 31*IntegerPow<2>(z)*(-30 + 198*IntegerPow<2>(z) - 495*IntegerPow<4>(z) + 407*IntegerPow<6>(z)))))/131072.;
    if constexpr(I == 390) return (-3*sqrt(98025655/(2.*pi))*z*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(256*IntegerPow<8>(y) + 512*IntegerPow<6>(y)*(-1 + IntegerPow<2>(z)) + 304*IntegerPow<4>(y)*IntegerPow<2>(-1 + IntegerPow<2>(z)) + 48*IntegerPow<2>(y)*IntegerPow<3>(-1 + IntegerPow<2>(z)) + IntegerPow<4>(-1 + IntegerPow<2>(z)))*(45 + 31*IntegerPow<2>(z)*(-60 + 594*IntegerPow<2>(z) - 1980*IntegerPow<4>(z) + 2035*IntegerPow<6>(z))))/65536.;
    if constexpr(I == 391) return (15*sqrt(58815393/pi)*(IntegerPow<11>(x) - 55*IntegerPow<9>(x)*IntegerPow<2>(y) + 330*IntegerPow<7>(x)*IntegerPow<4>(y) - 462*IntegerPow<5>(x)*IntegerPow<6>(y) + 165*IntegerPow<3>(x)*IntegerPow<8>(y) - 11*x*IntegerPow<10>(y))*(1 + 31*IntegerPow<2>(z)*(-4 + 66*IntegerPow<2>(z) - 308*IntegerPow<4>(z) + 407*IntegerPow<6>(z))))/131072.;
    if constexpr(I == 392) return (15*sqrt(1823277183/(2.*pi))*(IntegerPow<12>(x) - 66*IntegerPow<10>(x)*IntegerPow<2>(y) + 495*IntegerPow<8>(x)*IntegerPow<4>(y) - 924*IntegerPow<6>(x)*IntegerPow<6>(y) + 495*IntegerPow<4>(x)*IntegerPow<8>(y) - 66*IntegerPow<2>(x)*IntegerPow<10>(y) + IntegerPow<12>(y))*z*(-1 + 33*IntegerPow<2>(z) - 231*IntegerPow<4>(z) + 407*IntegerPow<6>(z)))/32768.;
    if constexpr(I == 393) return (15*sqrt(260468169/pi)*(IntegerPow<13>(x) - 78*IntegerPow<11>(x)*IntegerPow<2>(y) + 715*IntegerPow<9>(x)*IntegerPow<4>(y) - 1716*IntegerPow<7>(x)*IntegerPow<6>(y) + 1287*IntegerPow<5>(x)*IntegerPow<8>(y) - 286*IntegerPow<3>(x)*IntegerPow<10>(y) + 13*x*IntegerPow<12>(y))*(-1 + 99*IntegerPow<2>(z) - 1155*IntegerPow<4>(z) + 2849*IntegerPow<6>(z)))/262144.;
    if constexpr(I == 394) return (15*sqrt(2865149859/(2.*pi))*(IntegerPow<14>(x) - 91*IntegerPow<12>(x)*IntegerPow<2>(y) + 1001*IntegerPow<10>(x)*IntegerPow<4>(y) - 3003*IntegerPow<8>(x)*IntegerPow<6>(y) + 3003*IntegerPow<6>(x)*IntegerPow<8>(y) - 1001*IntegerPow<4>(x)*IntegerPow<10>(y) + 91*IntegerPow<2>(x)*IntegerPow<12>(y) - IntegerPow<14>(y))*z*(3 - 70*IntegerPow<2>(z) + 259*IntegerPow<4>(z)))/131072.;
    if constexpr(I == 395) return (3*sqrt(842691135/pi)*x*(IntegerPow<2>(x) - 3*IntegerPow<2>(y))*(IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + 5*IntegerPow<4>(y))*(IntegerPow<8>(x) - 92*IntegerPow<6>(x)*IntegerPow<2>(y) + 134*IntegerPow<4>(x)*IntegerPow<4>(y) - 28*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*(3 + 35*IntegerPow<2>(z)*(-6 + 37*IntegerPow<2>(z))))/262144.;
    if constexpr(I == 396) return (15*sqrt(1179767589/pi)*(IntegerPow<16>(x) - 120*IntegerPow<14>(x)*IntegerPow<2>(y) + 1820*IntegerPow<12>(x)*IntegerPow<4>(y) - 8008*IntegerPow<10>(x)*IntegerPow<6>(y) + 12870*IntegerPow<8>(x)*IntegerPow<8>(y) - 8008*IntegerPow<6>(x)*IntegerPow<10>(y) + 1820*IntegerPow<4>(x)*IntegerPow<12>(y) - 120*IntegerPow<2>(x)*IntegerPow<14>(y) + IntegerPow<16>(y))*z*(-3 + 37*IntegerPow<2>(z)))/131072.;
    if constexpr(I == 397) return (15*sqrt(393255863/pi)*(IntegerPow<17>(x) - 136*IntegerPow<15>(x)*IntegerPow<2>(y) + 2380*IntegerPow<13>(x)*IntegerPow<4>(y) - 12376*IntegerPow<11>(x)*IntegerPow<6>(y) + 24310*IntegerPow<9>(x)*IntegerPow<8>(y) - 19448*IntegerPow<7>(x)*IntegerPow<10>(y) + 6188*IntegerPow<5>(x)*IntegerPow<12>(y) - 680*IntegerPow<3>(x)*IntegerPow<14>(y) + 17*x*IntegerPow<16>(y))*(-1 + 37*IntegerPow<2>(z)))/262144.;
    if constexpr(I == 398) return (15*sqrt(14550466931/(2.*pi))*(IntegerPow<18>(x) - 153*IntegerPow<16>(x)*IntegerPow<2>(y) + 3060*IntegerPow<14>(x)*IntegerPow<4>(y) - 18564*IntegerPow<12>(x)*IntegerPow<6>(y) + 43758*IntegerPow<10>(x)*IntegerPow<8>(y) - 43758*IntegerPow<8>(x)*IntegerPow<10>(y) + 18564*IntegerPow<6>(x)*IntegerPow<12>(y) - 3060*IntegerPow<4>(x)*IntegerPow<14>(y) + 153*IntegerPow<2>(x)*IntegerPow<16>(y) - IntegerPow<18>(y))*z)/131072.;
    if constexpr(I == 399) return (15*sqrt(765814049/pi)*(IntegerPow<19>(x) - 171*IntegerPow<17>(x)*IntegerPow<2>(y) + 3876*IntegerPow<15>(x)*IntegerPow<4>(y) - 27132*IntegerPow<13>(x)*IntegerPow<6>(y) + 75582*IntegerPow<11>(x)*IntegerPow<8>(y) - 92378*IntegerPow<9>(x)*IntegerPow<10>(y) + 50388*IntegerPow<7>(x)*IntegerPow<12>(y) - 11628*IntegerPow<5>(x)*IntegerPow<14>(y) + 969*IntegerPow<3>(x)*IntegerPow<16>(y) - 19*x*IntegerPow<18>(y)))/262144.;
    if constexpr(I == 400) return (3*sqrt(156991880045/(2.*pi))*x*(x - y)*y*(x + y)*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(IntegerPow<4>(x) - 4*IntegerPow<3>(x)*y - 14*IntegerPow<2>(x)*IntegerPow<2>(y) - 4*x*IntegerPow<3>(y) + IntegerPow<4>(y))*(IntegerPow<4>(x) + 4*IntegerPow<3>(x)*y - 14*IntegerPow<2>(x)*IntegerPow<2>(y) + 4*x*IntegerPow<3>(y) + IntegerPow<4>(y))*(IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + 5*IntegerPow<4>(y)))/131072.;
    if constexpr(I == 401) return (-15*sqrt(31398376009/pi)*y*(-19*IntegerPow<18>(x) + 969*IntegerPow<16>(x)*IntegerPow<2>(y) - 11628*IntegerPow<14>(x)*IntegerPow<4>(y) + 50388*IntegerPow<12>(x)*IntegerPow<6>(y) - 92378*IntegerPow<10>(x)*IntegerPow<8>(y) + 75582*IntegerPow<8>(x)*IntegerPow<10>(y) - 27132*IntegerPow<6>(x)*IntegerPow<12>(y) + 3876*IntegerPow<4>(x)*IntegerPow<14>(y) - 171*IntegerPow<2>(x)*IntegerPow<16>(y) + IntegerPow<18>(y))*z)/262144.;
    if constexpr(I == 402) return (5*sqrt(7245779079/(2.*pi))*x*y*(IntegerPow<2>(x) - 3*IntegerPow<2>(y))*(3*IntegerPow<2>(x) - IntegerPow<2>(y))*(IntegerPow<6>(x) - 33*IntegerPow<4>(x)*IntegerPow<2>(y) + 27*IntegerPow<2>(x)*IntegerPow<4>(y) - 3*IntegerPow<6>(y))*(3*IntegerPow<6>(x) - 27*IntegerPow<4>(x)*IntegerPow<2>(y) + 33*IntegerPow<2>(x)*IntegerPow<4>(y) - IntegerPow<6>(y))*(-1 + 39*IntegerPow<2>(z)))/131072.;
    if constexpr(I == 403) return (15*sqrt(45889934167/pi)*y*(17*IntegerPow<16>(x) - 680*IntegerPow<14>(x)*IntegerPow<2>(y) + 6188*IntegerPow<12>(x)*IntegerPow<4>(y) - 19448*IntegerPow<10>(x)*IntegerPow<6>(y) + 24310*IntegerPow<8>(x)*IntegerPow<8>(y) - 12376*IntegerPow<6>(x)*IntegerPow<10>(y) + 2380*IntegerPow<4>(x)*IntegerPow<12>(y) - 136*IntegerPow<2>(x)*IntegerPow<14>(y) + IntegerPow<16>(y))*z*(-1 + 13*IntegerPow<2>(z)))/262144.;
    if constexpr(I == 404) return (-15*sqrt(1240268491/pi)*y*(-IntegerPow<15>(x) + 35*IntegerPow<13>(x)*IntegerPow<2>(y) - 273*IntegerPow<11>(x)*IntegerPow<4>(y) + 715*IntegerPow<9>(x)*IntegerPow<6>(y) - 715*IntegerPow<7>(x)*IntegerPow<8>(y) + 273*IntegerPow<5>(x)*IntegerPow<10>(y) - 35*IntegerPow<3>(x)*IntegerPow<12>(y) + x*IntegerPow<14>(y))*(1 - 74*IntegerPow<2>(z) + 481*IntegerPow<4>(z)))/32768.;
    if constexpr(I == 405) return (-3*sqrt(6201342455/pi)*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(IntegerPow<8>(x) - 28*IntegerPow<6>(x)*IntegerPow<2>(y) + 134*IntegerPow<4>(x)*IntegerPow<4>(y) - 92*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*z*(15 - 370*IntegerPow<2>(z) + 1443*IntegerPow<4>(z)))/262144.;
    if constexpr(I == 406) return (15*sqrt(531543639/(2.*pi))*x*y*(IntegerPow<6>(x) - 21*IntegerPow<4>(x)*IntegerPow<2>(y) + 35*IntegerPow<2>(x)*IntegerPow<4>(y) - 7*IntegerPow<6>(y))*(7*IntegerPow<6>(x) - 35*IntegerPow<4>(x)*IntegerPow<2>(y) + 21*IntegerPow<2>(x)*IntegerPow<4>(y) - IntegerPow<6>(y))*(-1 + 7*IntegerPow<2>(z)*(15 - 185*IntegerPow<2>(z) + 481*IntegerPow<4>(z))))/131072.;
    if constexpr(I == 407) return (15*sqrt(63253693041/pi)*y*(13*IntegerPow<12>(x) - 286*IntegerPow<10>(x)*IntegerPow<2>(y) + 1287*IntegerPow<8>(x)*IntegerPow<4>(y) - 1716*IntegerPow<6>(x)*IntegerPow<6>(y) + 715*IntegerPow<4>(x)*IntegerPow<8>(y) - 78*IntegerPow<2>(x)*IntegerPow<10>(y) + IntegerPow<12>(y))*z*(-1 + 35*IntegerPow<2>(z) - 259*IntegerPow<4>(z) + 481*IntegerPow<6>(z)))/262144.;
    if constexpr(I == 408) return (15*sqrt(1916778577/(2.*pi))*x*(x - y)*y*(x + y)*(IntegerPow<2>(x) - 3*IntegerPow<2>(y))*(3*IntegerPow<2>(x) - IntegerPow<2>(y))*(IntegerPow<2>(x) - 4*x*y + IntegerPow<2>(y))*(IntegerPow<2>(x) + 4*x*y + IntegerPow<2>(y))*(1 + 11*IntegerPow<2>(z)*(-12 + 210*IntegerPow<2>(z) - 1036*IntegerPow<4>(z) + 1443*IntegerPow<6>(z))))/131072.;
    if constexpr(I == 409) return (-15*sqrt(1916778577/pi)*y*(-11*IntegerPow<10>(x) + 165*IntegerPow<8>(x)*IntegerPow<2>(y) - 462*IntegerPow<6>(x)*IntegerPow<4>(y) + 330*IntegerPow<4>(x)*IntegerPow<6>(y) - 55*IntegerPow<2>(x)*IntegerPow<8>(y) + IntegerPow<10>(y))*z*(3 + 11*IntegerPow<2>(z)*(-12 + 126*IntegerPow<2>(z) - 444*IntegerPow<4>(z) + 481*IntegerPow<6>(z))))/131072.;
    if constexpr(I == 410) return (3*sqrt(309157835/(2.*pi))*x*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + 5*IntegerPow<4>(y))*(-3 + 465*IntegerPow<2>(z) + 341*IntegerPow<4>(z)*(-30 + 210*IntegerPow<2>(z) - 555*IntegerPow<4>(z) + 481*IntegerPow<6>(z))))/65536.;
    if constexpr(I == 411) return (5*sqrt(2040441711/pi)*y*(9*IntegerPow<8>(x) - 84*IntegerPow<6>(x)*IntegerPow<2>(y) + 126*IntegerPow<4>(x)*IntegerPow<4>(y) - 36*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*z*(-9 + 31*IntegerPow<2>(z)*(15 - 198*IntegerPow<2>(z) + 990*IntegerPow<4>(z) - 2035*IntegerPow<6>(z) + 1443*IntegerPow<8>(z))))/131072.;
    if constexpr(I == 412) return (-15*sqrt(23453353/pi)*y*(-IntegerPow<7>(x) + 7*IntegerPow<5>(x)*IntegerPow<2>(y) - 7*IntegerPow<3>(x)*IntegerPow<4>(y) + x*IntegerPow<6>(y))*(3 - 522*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(15 - 132*IntegerPow<2>(z) + 495*IntegerPow<4>(z) - 814*IntegerPow<6>(z) + 481*IntegerPow<8>(z))))/32768.;
    if constexpr(I == 413) return (-15*sqrt(43556227/pi)*y*(-7*IntegerPow<6>(x) + 35*IntegerPow<4>(x)*IntegerPow<2>(y) - 21*IntegerPow<2>(x)*IntegerPow<4>(y) + IntegerPow<6>(y))*z*(21 + 29*IntegerPow<2>(z)*(-42 + 31*IntegerPow<2>(z)*(21 - 132*IntegerPow<2>(z) + 7*IntegerPow<4>(z)*(55 + 37*IntegerPow<2>(z)*(-2 + IntegerPow<2>(z)))))))/131072.;
    if constexpr(I == 414) return (5*sqrt(914680767/(2.*pi))*y*(3*IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 3*x*IntegerPow<4>(y))*(-1 + 189*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(-63 + 31*IntegerPow<2>(z)*(21 - 99*IntegerPow<2>(z) + 231*IntegerPow<4>(z) - 259*IntegerPow<6>(z) + 111*IntegerPow<8>(z)))))/65536.;
    if constexpr(I == 415) return (3*sqrt(117266765/pi)*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*z*(-65 + 4095*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(-819 + 31*IntegerPow<2>(z)*(195 - 715*IntegerPow<2>(z) + 1365*IntegerPow<4>(z) - 1295*IntegerPow<6>(z) + 481*IntegerPow<8>(z)))))/131072.;
    if constexpr(I == 416) return (-3*sqrt(117266765/pi)*x*y*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(13 + 5*IntegerPow<2>(z)*(-520 + 3*IntegerPow<2>(z)*(5460 - 63336*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(390 - 1144*IntegerPow<2>(z) + 1820*IntegerPow<4>(z) - 1480*IntegerPow<6>(z) + 481*IntegerPow<8>(z))))))/131072.;
    if constexpr(I == 417) return -(sqrt(20694135/(2.*pi))*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*z*(663 + 5*IntegerPow<2>(z)*(-8840 + 3*IntegerPow<2>(z)*(55692 - 461448*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(2210 - 5304*IntegerPow<2>(z) + 7140*IntegerPow<4>(z) - 5032*IntegerPow<6>(z) + 1443*IntegerPow<8>(z))))))/131072.;
    if constexpr(I == 418) return (sqrt(899745/pi)*x*y*(-221 + 69*IntegerPow<2>(z)*(663 + 5*IntegerPow<2>(z)*(-4420 + 55692*IntegerPow<2>(z) - 346086*IntegerPow<4>(z) + 899*IntegerPow<6>(z)*(1326 - 2652*IntegerPow<2>(z) + 3060*IntegerPow<4>(z) - 1887*IntegerPow<6>(z) + 481*IntegerPow<8>(z))))))/131072.;
    if constexpr(I == 419) return (sqrt(4305/(2.*pi))*y*z*(-46189 + 69*IntegerPow<2>(z)*(46189 - 923780*IntegerPow<2>(z) + 8314020*IntegerPow<4>(z) + 145*IntegerPow<6>(z)*(-277134 + 781014*IntegerPow<2>(z) + 341*IntegerPow<4>(z)*(-3876 + 3876*IntegerPow<2>(z) - 2109*IntegerPow<4>(z) + 481*IntegerPow<6>(z))))))/131072.;
    if constexpr(I == 420) return (sqrt(41/pi)*(46189 - 9699690*IntegerPow<2>(z) + 345*IntegerPow<4>(z)*(969969 - 12932920*IntegerPow<2>(z) + 87297210*IntegerPow<4>(z) - 337549212*IntegerPow<6>(z) + 899*IntegerPow<8>(z)*(881790 - 1279080*IntegerPow<2>(z) + 77*IntegerPow<4>(z)*(14535 - 7030*IntegerPow<2>(z) + 1443*IntegerPow<4>(z))))))/524288.;
    if constexpr(I == 421) return (sqrt(4305/(2.*pi))*x*z*(-46189 + 69*IntegerPow<2>(z)*(46189 - 923780*IntegerPow<2>(z) + 8314020*IntegerPow<4>(z) + 145*IntegerPow<6>(z)*(-277134 + 781014*IntegerPow<2>(z) + 341*IntegerPow<4>(z)*(-3876 + 3876*IntegerPow<2>(z) - 2109*IntegerPow<4>(z) + 481*IntegerPow<6>(z))))))/131072.;
    if constexpr(I == 422) return -(sqrt(899745/pi)*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(-221 + 69*IntegerPow<2>(z)*(663 + 5*IntegerPow<2>(z)*(-4420 + 55692*IntegerPow<2>(z) - 346086*IntegerPow<4>(z) + 899*IntegerPow<6>(z)*(1326 - 2652*IntegerPow<2>(z) + 3060*IntegerPow<4>(z) - 1887*IntegerPow<6>(z) + 481*IntegerPow<8>(z))))))/262144.;
    if constexpr(I == 423) return (sqrt(20694135/(2.*pi))*(IntegerPow<3>(x) - 3*x*IntegerPow<2>(y))*z*(663 + 5*IntegerPow<2>(z)*(-8840 + 3*IntegerPow<2>(z)*(55692 - 461448*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(2210 - 5304*IntegerPow<2>(z) + 7140*IntegerPow<4>(z) - 5032*IntegerPow<6>(z) + 1443*IntegerPow<8>(z))))))/131072.;
    if constexpr(I == 424) return (3*sqrt(117266765/pi)*(8*IntegerPow<4>(y) + 8*IntegerPow<2>(y)*(-1 + IntegerPow<2>(z)) + IntegerPow<2>(-1 + IntegerPow<2>(z)))*(13 + 5*IntegerPow<2>(z)*(-520 + 3*IntegerPow<2>(z)*(5460 - 63336*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(390 - 1144*IntegerPow<2>(z) + 1820*IntegerPow<4>(z) - 1480*IntegerPow<6>(z) + 481*IntegerPow<8>(z))))))/524288.;
    if constexpr(I == 425) return (3*sqrt(117266765/pi)*(IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 5*x*IntegerPow<4>(y))*z*(-65 + 4095*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(-819 + 31*IntegerPow<2>(z)*(195 - 715*IntegerPow<2>(z) + 1365*IntegerPow<4>(z) - 1295*IntegerPow<6>(z) + 481*IntegerPow<8>(z)))))/131072.;
    if constexpr(I == 426) return (5*sqrt(914680767/(2.*pi))*(IntegerPow<6>(x) - 15*IntegerPow<4>(x)*IntegerPow<2>(y) + 15*IntegerPow<2>(x)*IntegerPow<4>(y) - IntegerPow<6>(y))*(-1 + 189*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(-63 + 31*IntegerPow<2>(z)*(21 - 99*IntegerPow<2>(z) + 231*IntegerPow<4>(z) - 259*IntegerPow<6>(z) + 111*IntegerPow<8>(z)))))/131072.;
    if constexpr(I == 427) return (15*sqrt(43556227/pi)*(IntegerPow<7>(x) - 21*IntegerPow<5>(x)*IntegerPow<2>(y) + 35*IntegerPow<3>(x)*IntegerPow<4>(y) - 7*x*IntegerPow<6>(y))*z*(21 + 29*IntegerPow<2>(z)*(-42 + 31*IntegerPow<2>(z)*(21 - 132*IntegerPow<2>(z) + 7*IntegerPow<4>(z)*(55 + 37*IntegerPow<2>(z)*(-2 + IntegerPow<2>(z)))))))/131072.;
    if constexpr(I == 428) return (15*sqrt(23453353/pi)*(128*IntegerPow<8>(y) + 256*IntegerPow<6>(y)*(-1 + IntegerPow<2>(z)) + 160*IntegerPow<4>(y)*IntegerPow<2>(-1 + IntegerPow<2>(z)) + 32*IntegerPow<2>(y)*IntegerPow<3>(-1 + IntegerPow<2>(z)) + IntegerPow<4>(-1 + IntegerPow<2>(z)))*(3 - 522*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(15 - 132*IntegerPow<2>(z) + 495*IntegerPow<4>(z) - 814*IntegerPow<6>(z) + 481*IntegerPow<8>(z))))/262144.;
    if constexpr(I == 429) return (5*sqrt(2040441711/pi)*(IntegerPow<9>(x) - 36*IntegerPow<7>(x)*IntegerPow<2>(y) + 126*IntegerPow<5>(x)*IntegerPow<4>(y) - 84*IntegerPow<3>(x)*IntegerPow<6>(y) + 9*x*IntegerPow<8>(y))*z*(-9 + 31*IntegerPow<2>(z)*(15 - 198*IntegerPow<2>(z) + 990*IntegerPow<4>(z) - 2035*IntegerPow<6>(z) + 1443*IntegerPow<8>(z))))/131072.;
    if constexpr(I == 430) return (-3*sqrt(309157835/(2.*pi))*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(256*IntegerPow<8>(y) + 512*IntegerPow<6>(y)*(-1 + IntegerPow<2>(z)) + 304*IntegerPow<4>(y)*IntegerPow<2>(-1 + IntegerPow<2>(z)) + 48*IntegerPow<2>(y)*IntegerPow<3>(-1 + IntegerPow<2>(z)) + IntegerPow<4>(-1 + IntegerPow<2>(z)))*(-3 + 465*IntegerPow<2>(z) + 341*IntegerPow<4>(z)*(-30 + 210*IntegerPow<2>(z) - 555*IntegerPow<4>(z) + 481*IntegerPow<6>(z))))/131072.;
    if constexpr(I == 431) return (15*sqrt(1916778577/pi)*(IntegerPow<11>(x) - 55*IntegerPow<9>(x)*IntegerPow<2>(y) + 330*IntegerPow<7>(x)*IntegerPow<4>(y) - 462*IntegerPow<5>(x)*IntegerPow<6>(y) + 165*IntegerPow<3>(x)*IntegerPow<8>(y) - 11*x*IntegerPow<10>(y))*z*(3 + 11*IntegerPow<2>(z)*(-12 + 126*IntegerPow<2>(z) - 444*IntegerPow<4>(z) + 481*IntegerPow<6>(z))))/131072.;
    if constexpr(I == 432) return (15*sqrt(1916778577/(2.*pi))*(2048*IntegerPow<12>(y) + 6144*IntegerPow<10>(y)*(-1 + IntegerPow<2>(z)) + 6912*IntegerPow<8>(y)*IntegerPow<2>(-1 + IntegerPow<2>(z)) + 3584*IntegerPow<6>(y)*IntegerPow<3>(-1 + IntegerPow<2>(z)) + 840*IntegerPow<4>(y)*IntegerPow<4>(-1 + IntegerPow<2>(z)) + 72*IntegerPow<2>(y)*IntegerPow<5>(-1 + IntegerPow<2>(z)) + IntegerPow<6>(-1 + IntegerPow<2>(z)))*(1 + 11*IntegerPow<2>(z)*(-12 + 210*IntegerPow<2>(z) - 1036*IntegerPow<4>(z) + 1443*IntegerPow<6>(z))))/524288.;
    if constexpr(I == 433) return (15*sqrt(63253693041/pi)*(IntegerPow<13>(x) - 78*IntegerPow<11>(x)*IntegerPow<2>(y) + 715*IntegerPow<9>(x)*IntegerPow<4>(y) - 1716*IntegerPow<7>(x)*IntegerPow<6>(y) + 1287*IntegerPow<5>(x)*IntegerPow<8>(y) - 286*IntegerPow<3>(x)*IntegerPow<10>(y) + 13*x*IntegerPow<12>(y))*z*(-1 + 35*IntegerPow<2>(z) - 259*IntegerPow<4>(z) + 481*IntegerPow<6>(z)))/262144.;
    if constexpr(I == 434) return (15*sqrt(531543639/(2.*pi))*(IntegerPow<14>(x) - 91*IntegerPow<12>(x)*IntegerPow<2>(y) + 1001*IntegerPow<10>(x)*IntegerPow<4>(y) - 3003*IntegerPow<8>(x)*IntegerPow<6>(y) + 3003*IntegerPow<6>(x)*IntegerPow<8>(y) - 1001*IntegerPow<4>(x)*IntegerPow<10>(y) + 91*IntegerPow<2>(x)*IntegerPow<12>(y) - IntegerPow<14>(y))*(-1 + 7*IntegerPow<2>(z)*(15 - 185*IntegerPow<2>(z) + 481*IntegerPow<4>(z))))/262144.;
    if constexpr(I == 435) return (3*sqrt(6201342455/pi)*x*(IntegerPow<2>(x) - 3*IntegerPow<2>(y))*(IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + 5*IntegerPow<4>(y))*(IntegerPow<8>(x) - 92*IntegerPow<6>(x)*IntegerPow<2>(y) + 134*IntegerPow<4>(x)*IntegerPow<4>(y) - 28*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*z*(15 - 370*IntegerPow<2>(z) + 1443*IntegerPow<4>(z)))/262144.;
    if constexpr(I == 436) return (15*sqrt(1240268491/pi)*(IntegerPow<16>(x) - 120*IntegerPow<14>(x)*IntegerPow<2>(y) + 1820*IntegerPow<12>(x)*IntegerPow<4>(y) - 8008*IntegerPow<10>(x)*IntegerPow<6>(y) + 12870*IntegerPow<8>(x)*IntegerPow<8>(y) - 8008*IntegerPow<6>(x)*IntegerPow<10>(y) + 1820*IntegerPow<4>(x)*IntegerPow<12>(y) - 120*IntegerPow<2>(x)*IntegerPow<14>(y) + IntegerPow<16>(y))*(1 - 74*IntegerPow<2>(z) + 481*IntegerPow<4>(z)))/524288.;
    if constexpr(I == 437) return (15*sqrt(45889934167/pi)*(IntegerPow<17>(x) - 136*IntegerPow<15>(x)*IntegerPow<2>(y) + 2380*IntegerPow<13>(x)*IntegerPow<4>(y) - 12376*IntegerPow<11>(x)*IntegerPow<6>(y) + 24310*IntegerPow<9>(x)*IntegerPow<8>(y) - 19448*IntegerPow<7>(x)*IntegerPow<10>(y) + 6188*IntegerPow<5>(x)*IntegerPow<12>(y) - 680*IntegerPow<3>(x)*IntegerPow<14>(y) + 17*x*IntegerPow<16>(y))*z*(-1 + 13*IntegerPow<2>(z)))/262144.;
    if constexpr(I == 438) return (5*sqrt(7245779079/(2.*pi))*(IntegerPow<18>(x) - 153*IntegerPow<16>(x)*IntegerPow<2>(y) + 3060*IntegerPow<14>(x)*IntegerPow<4>(y) - 18564*IntegerPow<12>(x)*IntegerPow<6>(y) + 43758*IntegerPow<10>(x)*IntegerPow<8>(y) - 43758*IntegerPow<8>(x)*IntegerPow<10>(y) + 18564*IntegerPow<6>(x)*IntegerPow<12>(y) - 3060*IntegerPow<4>(x)*IntegerPow<14>(y) + 153*IntegerPow<2>(x)*IntegerPow<16>(y) - IntegerPow<18>(y))*(-1 + 39*IntegerPow<2>(z)))/262144.;
    if constexpr(I == 439) return (15*sqrt(31398376009/pi)*(IntegerPow<19>(x) - 171*IntegerPow<17>(x)*IntegerPow<2>(y) + 3876*IntegerPow<15>(x)*IntegerPow<4>(y) - 27132*IntegerPow<13>(x)*IntegerPow<6>(y) + 75582*IntegerPow<11>(x)*IntegerPow<8>(y) - 92378*IntegerPow<9>(x)*IntegerPow<10>(y) + 50388*IntegerPow<7>(x)*IntegerPow<12>(y) - 11628*IntegerPow<5>(x)*IntegerPow<14>(y) + 969*IntegerPow<3>(x)*IntegerPow<16>(y) - 19*x*IntegerPow<18>(y))*z)/262144.;
    if constexpr(I == 440) return (3*sqrt(156991880045/(2.*pi))*(IntegerPow<20>(x) - 190*IntegerPow<18>(x)*IntegerPow<2>(y) + 4845*IntegerPow<16>(x)*IntegerPow<4>(y) - 38760*IntegerPow<14>(x)*IntegerPow<6>(y) + 125970*IntegerPow<12>(x)*IntegerPow<8>(y) - 184756*IntegerPow<10>(x)*IntegerPow<10>(y) + 125970*IntegerPow<8>(x)*IntegerPow<12>(y) - 38760*IntegerPow<6>(x)*IntegerPow<14>(y) + 4845*IntegerPow<4>(x)*IntegerPow<16>(y) - 190*IntegerPow<2>(x)*IntegerPow<18>(y) + IntegerPow<20>(y)))/524288.;
    if constexpr(I == 441) return (sqrt(2893136075115/pi)*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*(-7*IntegerPow<6>(x) + 35*IntegerPow<4>(x)*IntegerPow<2>(y) - 21*IntegerPow<2>(x)*IntegerPow<4>(y) + IntegerPow<6>(y))*(IntegerPow<12>(x) - 58*IntegerPow<10>(x)*IntegerPow<2>(y) + 655*IntegerPow<8>(x)*IntegerPow<4>(y) - 1772*IntegerPow<6>(x)*IntegerPow<6>(y) + 1423*IntegerPow<4>(x)*IntegerPow<8>(y) - 186*IntegerPow<2>(x)*IntegerPow<10>(y) + IntegerPow<12>(y)))/1.048576e6;
    if constexpr(I == 442) return (3*sqrt(6750650841935/(2.*pi))*x*(x - y)*y*(x + y)*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(IntegerPow<4>(x) - 4*IntegerPow<3>(x)*y - 14*IntegerPow<2>(x)*IntegerPow<2>(y) - 4*x*IntegerPow<3>(y) + IntegerPow<4>(y))*(IntegerPow<4>(x) + 4*IntegerPow<3>(x)*y - 14*IntegerPow<2>(x)*IntegerPow<2>(y) + 4*x*IntegerPow<3>(y) + IntegerPow<4>(y))*(IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + 5*IntegerPow<4>(y))*z)/131072.;
    if constexpr(I == 443) return (-3*sqrt(164650020535/pi)*y*(-19*IntegerPow<18>(x) + 969*IntegerPow<16>(x)*IntegerPow<2>(y) - 11628*IntegerPow<14>(x)*IntegerPow<4>(y) + 50388*IntegerPow<12>(x)*IntegerPow<6>(y) - 92378*IntegerPow<10>(x)*IntegerPow<8>(y) + 75582*IntegerPow<8>(x)*IntegerPow<10>(y) - 27132*IntegerPow<6>(x)*IntegerPow<12>(y) + 3876*IntegerPow<4>(x)*IntegerPow<14>(y) - 171*IntegerPow<2>(x)*IntegerPow<16>(y) + IntegerPow<18>(y))*(-1 + 41*IntegerPow<2>(z)))/1.048576e6;
    if constexpr(I == 444) return (5*sqrt(98790012321/(2.*pi))*x*y*(IntegerPow<2>(x) - 3*IntegerPow<2>(y))*(3*IntegerPow<2>(x) - IntegerPow<2>(y))*(IntegerPow<6>(x) - 33*IntegerPow<4>(x)*IntegerPow<2>(y) + 27*IntegerPow<2>(x)*IntegerPow<4>(y) - 3*IntegerPow<6>(y))*(3*IntegerPow<6>(x) - 27*IntegerPow<4>(x)*IntegerPow<2>(y) + 33*IntegerPow<2>(x)*IntegerPow<4>(y) - IntegerPow<6>(y))*z*(-3 + 41*IntegerPow<2>(z)))/131072.;
    if constexpr(I == 445) return (15*sqrt(2533077239/(2.*pi))*y*(17*IntegerPow<16>(x) - 680*IntegerPow<14>(x)*IntegerPow<2>(y) + 6188*IntegerPow<12>(x)*IntegerPow<4>(y) - 19448*IntegerPow<10>(x)*IntegerPow<6>(y) + 24310*IntegerPow<8>(x)*IntegerPow<8>(y) - 12376*IntegerPow<6>(x)*IntegerPow<10>(y) + 2380*IntegerPow<4>(x)*IntegerPow<12>(y) - 136*IntegerPow<2>(x)*IntegerPow<14>(y) + IntegerPow<16>(y))*(1 - 78*IntegerPow<2>(z) + 533*IntegerPow<4>(z)))/524288.;
    if constexpr(I == 446) return (-3*sqrt(240642337705/pi)*y*(-IntegerPow<15>(x) + 35*IntegerPow<13>(x)*IntegerPow<2>(y) - 273*IntegerPow<11>(x)*IntegerPow<4>(y) + 715*IntegerPow<9>(x)*IntegerPow<6>(y) - 715*IntegerPow<7>(x)*IntegerPow<8>(y) + 273*IntegerPow<5>(x)*IntegerPow<10>(y) - 35*IntegerPow<3>(x)*IntegerPow<12>(y) + x*IntegerPow<14>(y))*z*(5 - 130*IntegerPow<2>(z) + 533*IntegerPow<4>(z)))/32768.;
    if constexpr(I == 447) return -(sqrt(19511540895/(2.*pi))*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(IntegerPow<8>(x) - 28*IntegerPow<6>(x)*IntegerPow<2>(y) + 134*IntegerPow<4>(x)*IntegerPow<4>(y) - 92*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*(-5 + 37*IntegerPow<2>(z)*(15 - 195*IntegerPow<2>(z) + 533*IntegerPow<4>(z))))/524288.;
    if constexpr(I == 448) return (3*sqrt(2787362985/(2.*pi))*x*y*(IntegerPow<6>(x) - 21*IntegerPow<4>(x)*IntegerPow<2>(y) + 35*IntegerPow<2>(x)*IntegerPow<4>(y) - 7*IntegerPow<6>(y))*(7*IntegerPow<6>(x) - 35*IntegerPow<4>(x)*IntegerPow<2>(y) + 21*IntegerPow<2>(x)*IntegerPow<4>(y) - IntegerPow<6>(y))*z*(-35 + 37*IntegerPow<2>(z)*(35 - 273*IntegerPow<2>(z) + 533*IntegerPow<4>(z))))/131072.;
    if constexpr(I == 449) return (15*sqrt(3902308179/pi)*y*(13*IntegerPow<12>(x) - 286*IntegerPow<10>(x)*IntegerPow<2>(y) + 1287*IntegerPow<8>(x)*IntegerPow<4>(y) - 1716*IntegerPow<6>(x)*IntegerPow<6>(y) + 715*IntegerPow<4>(x)*IntegerPow<8>(y) - 78*IntegerPow<2>(x)*IntegerPow<10>(y) + IntegerPow<12>(y))*(1 - 140*IntegerPow<2>(z) + 37*IntegerPow<4>(z)*(70 - 364*IntegerPow<2>(z) + 533*IntegerPow<4>(z))))/1.048576e6;
    if constexpr(I == 450) return (5*sqrt(66339239043/(2.*pi))*x*(x - y)*y*(x + y)*(IntegerPow<2>(x) - 3*IntegerPow<2>(y))*(3*IntegerPow<2>(x) - IntegerPow<2>(y))*(IntegerPow<2>(x) - 4*x*y + IntegerPow<2>(y))*(IntegerPow<2>(x) + 4*x*y + IntegerPow<2>(y))*z*(9 - 420*IntegerPow<2>(z) + 37*IntegerPow<4>(z)*(126 - 468*IntegerPow<2>(z) + 533*IntegerPow<4>(z))))/131072.;
    if constexpr(I == 451) return (-3*sqrt(10051399855/pi)*y*(-11*IntegerPow<10>(x) + 165*IntegerPow<8>(x)*IntegerPow<2>(y) - 462*IntegerPow<6>(x)*IntegerPow<4>(y) + 330*IntegerPow<4>(x)*IntegerPow<6>(y) - 55*IntegerPow<2>(x)*IntegerPow<8>(y) + IntegerPow<10>(y))*(-3 + 11*IntegerPow<2>(z)*(45 - 1050*IntegerPow<2>(z) + 37*IntegerPow<4>(z)*(210 - 585*IntegerPow<2>(z) + 533*IntegerPow<4>(z)))))/1.048576e6;
    if constexpr(I == 452) return (3*sqrt(110565398405/(2.*pi))*x*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + 5*IntegerPow<4>(y))*z*(-3 + 165*IntegerPow<2>(z) - 2310*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(330 - 715*IntegerPow<2>(z) + 533*IntegerPow<4>(z))))/65536.;
    if constexpr(I == 453) return (sqrt(10699877265/(2.*pi))*y*(9*IntegerPow<8>(x) - 84*IntegerPow<6>(x)*IntegerPow<2>(y) + 126*IntegerPow<4>(x)*IntegerPow<4>(y) - 36*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*(3 + 31*IntegerPow<2>(z)*(-18 + 495*IntegerPow<2>(z) - 4620*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(495 - 858*IntegerPow<2>(z) + 533*IntegerPow<4>(z)))))/262144.;
    if constexpr(I == 454) return (-15*sqrt(9273226963/pi)*y*(-IntegerPow<7>(x) + 7*IntegerPow<5>(x)*IntegerPow<2>(y) - 7*IntegerPow<3>(x)*IntegerPow<4>(y) + x*IntegerPow<6>(y))*z*(3 + 31*IntegerPow<2>(z)*(-6 + 99*IntegerPow<2>(z) - 660*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(55 - 78*IntegerPow<2>(z) + 41*IntegerPow<4>(z)))))/32768.;
    if constexpr(I == 455) return (-15*sqrt(45680921/(2.*pi))*y*(-7*IntegerPow<6>(x) + 35*IntegerPow<4>(x)*IntegerPow<2>(y) - 21*IntegerPow<2>(x)*IntegerPow<4>(y) + IntegerPow<6>(y))*(-3 + 609*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(-21 + 231*IntegerPow<2>(z) - 1155*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(77 - 91*IntegerPow<2>(z) + 41*IntegerPow<4>(z)))))/262144.;
    if constexpr(I == 456) return (sqrt(4796496705/(2.*pi))*y*(3*IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 3*x*IntegerPow<4>(y))*z*(-45 + 3045*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(-63 + 495*IntegerPow<2>(z) - 1925*IntegerPow<4>(z) + 3885*IntegerPow<6>(z) - 3885*IntegerPow<8>(z) + 1517*IntegerPow<10>(z))))/65536.;
    if constexpr(I == 457) return (3*sqrt(1598832235/(2.*pi))*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(5 - 1080*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(420 + 31*IntegerPow<2>(z)*(-168 + 990*IntegerPow<2>(z) - 3080*IntegerPow<4>(z) + 5180*IntegerPow<6>(z) - 4440*IntegerPow<8>(z) + 1517*IntegerPow<10>(z)))))/524288.;
    if constexpr(I == 458) return (-3*sqrt(7234535/pi)*x*y*z*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(1105 + 3*IntegerPow<2>(z)*(-26520 + 538356*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(-5304 + 24310*IntegerPow<2>(z) - 61880*IntegerPow<4>(z) + 88060*IntegerPow<6>(z) - 65416*IntegerPow<8>(z) + 19721*IntegerPow<10>(z)))))/131072.;
    if constexpr(I == 459) return -(sqrt(7234535/(2.*pi))*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*(-221 + 15*IntegerPow<2>(z)*(3315 - 119340*IntegerPow<2>(z) + 29*IntegerPow<4>(z)*(55692 + 31*IntegerPow<2>(z)*(-11934 + 43758*IntegerPow<2>(z) - 92820*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(3060 - 1989*IntegerPow<2>(z) + 533*IntegerPow<4>(z)))))))/524288.;
    if constexpr(I == 460) return (sqrt(1142295/pi)*x*y*z*(-4199 + 15*IntegerPow<2>(z)*(20995 - 453492*IntegerPow<2>(z) + 4383756*IntegerPow<4>(z) + 899*IntegerPow<6>(z)*(-25194 + 75582*IntegerPow<2>(z) - 135660*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(3876 - 2223*IntegerPow<2>(z) + 533*IntegerPow<4>(z))))))/131072.;
    if constexpr(I == 461) return (sqrt(9933/pi)*y*(4199 + 115*IntegerPow<2>(z)*(-8398 + 314925*IntegerPow<2>(z) - 4534920*IntegerPow<4>(z) + 87*IntegerPow<6>(z)*(377910 + 31*IntegerPow<2>(z)*(-50388 + 125970*IntegerPow<2>(z) - 193800*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(4845 - 2470*IntegerPow<2>(z) + 533*IntegerPow<4>(z)))))))/524288.;
    if constexpr(I == 462) return (sqrt(43/pi)*z*(969969 + 115*IntegerPow<2>(z)*(-646646 + 3*IntegerPow<2>(z)*(4849845 - 49884120*IntegerPow<2>(z) + 29*IntegerPow<4>(z)*(9699690 - 32802588*IntegerPow<2>(z) + 341*IntegerPow<4>(z)*(203490 - 271320*IntegerPow<2>(z) + 37*IntegerPow<4>(z)*(5985 - 2730*IntegerPow<2>(z) + 533*IntegerPow<4>(z))))))))/524288.;
    if constexpr(I == 463) return (sqrt(9933/pi)*x*(4199 + 115*IntegerPow<2>(z)*(-8398 + 314925*IntegerPow<2>(z) - 4534920*IntegerPow<4>(z) + 87*IntegerPow<6>(z)*(377910 + 31*IntegerPow<2>(z)*(-50388 + 125970*IntegerPow<2>(z) - 193800*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(4845 - 2470*IntegerPow<2>(z) + 533*IntegerPow<4>(z)))))))/524288.;
    if constexpr(I == 464) return -(sqrt(1142295/pi)*z*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(-4199 + 15*IntegerPow<2>(z)*(20995 - 453492*IntegerPow<2>(z) + 4383756*IntegerPow<4>(z) + 899*IntegerPow<6>(z)*(-25194 + 75582*IntegerPow<2>(z) - 135660*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(3876 - 2223*IntegerPow<2>(z) + 533*IntegerPow<4>(z))))))/262144.;
    if constexpr(I == 465) return (sqrt(7234535/(2.*pi))*(IntegerPow<3>(x) - 3*x*IntegerPow<2>(y))*(-221 + 15*IntegerPow<2>(z)*(3315 - 119340*IntegerPow<2>(z) + 29*IntegerPow<4>(z)*(55692 + 31*IntegerPow<2>(z)*(-11934 + 43758*IntegerPow<2>(z) - 92820*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(3060 - 1989*IntegerPow<2>(z) + 533*IntegerPow<4>(z)))))))/524288.;
    if constexpr(I == 466) return (3*sqrt(7234535/pi)*z*(8*IntegerPow<4>(y) + 8*IntegerPow<2>(y)*(-1 + IntegerPow<2>(z)) + IntegerPow<2>(-1 + IntegerPow<2>(z)))*(1105 + 3*IntegerPow<2>(z)*(-26520 + 538356*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(-5304 + 24310*IntegerPow<2>(z) - 61880*IntegerPow<4>(z) + 88060*IntegerPow<6>(z) - 65416*IntegerPow<8>(z) + 19721*IntegerPow<10>(z)))))/524288.;
    if constexpr(I == 467) return (3*sqrt(1598832235/(2.*pi))*(IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 5*x*IntegerPow<4>(y))*(5 - 1080*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(420 + 31*IntegerPow<2>(z)*(-168 + 990*IntegerPow<2>(z) - 3080*IntegerPow<4>(z) + 5180*IntegerPow<6>(z) - 4440*IntegerPow<8>(z) + 1517*IntegerPow<10>(z)))))/524288.;
    if constexpr(I == 468) return (sqrt(4796496705/(2.*pi))*(IntegerPow<6>(x) - 15*IntegerPow<4>(x)*IntegerPow<2>(y) + 15*IntegerPow<2>(x)*IntegerPow<4>(y) - IntegerPow<6>(y))*z*(-45 + 3045*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(-63 + 495*IntegerPow<2>(z) - 1925*IntegerPow<4>(z) + 3885*IntegerPow<6>(z) - 3885*IntegerPow<8>(z) + 1517*IntegerPow<10>(z))))/131072.;
    if constexpr(I == 469) return (15*sqrt(45680921/(2.*pi))*(IntegerPow<7>(x) - 21*IntegerPow<5>(x)*IntegerPow<2>(y) + 35*IntegerPow<3>(x)*IntegerPow<4>(y) - 7*x*IntegerPow<6>(y))*(-3 + 609*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(-21 + 231*IntegerPow<2>(z) - 1155*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(77 - 91*IntegerPow<2>(z) + 41*IntegerPow<4>(z)))))/262144.;
    if constexpr(I == 470) return (15*sqrt(9273226963/pi)*z*(128*IntegerPow<8>(y) + 256*IntegerPow<6>(y)*(-1 + IntegerPow<2>(z)) + 160*IntegerPow<4>(y)*IntegerPow<2>(-1 + IntegerPow<2>(z)) + 32*IntegerPow<2>(y)*IntegerPow<3>(-1 + IntegerPow<2>(z)) + IntegerPow<4>(-1 + IntegerPow<2>(z)))*(3 + 31*IntegerPow<2>(z)*(-6 + 99*IntegerPow<2>(z) - 660*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(55 - 78*IntegerPow<2>(z) + 41*IntegerPow<4>(z)))))/262144.;
    if constexpr(I == 471) return (sqrt(10699877265/(2.*pi))*(IntegerPow<9>(x) - 36*IntegerPow<7>(x)*IntegerPow<2>(y) + 126*IntegerPow<5>(x)*IntegerPow<4>(y) - 84*IntegerPow<3>(x)*IntegerPow<6>(y) + 9*x*IntegerPow<8>(y))*(3 + 31*IntegerPow<2>(z)*(-18 + 495*IntegerPow<2>(z) - 4620*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(495 - 858*IntegerPow<2>(z) + 533*IntegerPow<4>(z)))))/262144.;
    if constexpr(I == 472) return (-3*sqrt(110565398405/(2.*pi))*z*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(256*IntegerPow<8>(y) + 512*IntegerPow<6>(y)*(-1 + IntegerPow<2>(z)) + 304*IntegerPow<4>(y)*IntegerPow<2>(-1 + IntegerPow<2>(z)) + 48*IntegerPow<2>(y)*IntegerPow<3>(-1 + IntegerPow<2>(z)) + IntegerPow<4>(-1 + IntegerPow<2>(z)))*(-3 + 165*IntegerPow<2>(z) - 2310*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(330 - 715*IntegerPow<2>(z) + 533*IntegerPow<4>(z))))/131072.;
    if constexpr(I == 473) return (3*sqrt(10051399855/pi)*(IntegerPow<11>(x) - 55*IntegerPow<9>(x)*IntegerPow<2>(y) + 330*IntegerPow<7>(x)*IntegerPow<4>(y) - 462*IntegerPow<5>(x)*IntegerPow<6>(y) + 165*IntegerPow<3>(x)*IntegerPow<8>(y) - 11*x*IntegerPow<10>(y))*(-3 + 11*IntegerPow<2>(z)*(45 - 1050*IntegerPow<2>(z) + 37*IntegerPow<4>(z)*(210 - 585*IntegerPow<2>(z) + 533*IntegerPow<4>(z)))))/1.048576e6;
    if constexpr(I == 474) return (5*sqrt(66339239043/(2.*pi))*z*(8*IntegerPow<4>(y) + 8*IntegerPow<2>(y)*(-1 + IntegerPow<2>(z)) + IntegerPow<2>(-1 + IntegerPow<2>(z)))*(256*IntegerPow<8>(y) + 512*IntegerPow<6>(y)*(-1 + IntegerPow<2>(z)) + 320*IntegerPow<4>(y)*IntegerPow<2>(-1 + IntegerPow<2>(z)) + 64*IntegerPow<2>(y)*IntegerPow<3>(-1 + IntegerPow<2>(z)) + IntegerPow<4>(-1 + IntegerPow<2>(z)))*(9 - 420*IntegerPow<2>(z) + 37*IntegerPow<4>(z)*(126 - 468*IntegerPow<2>(z) + 533*IntegerPow<4>(z))))/524288.;
    if constexpr(I == 475) return (15*sqrt(3902308179/pi)*(IntegerPow<13>(x) - 78*IntegerPow<11>(x)*IntegerPow<2>(y) + 715*IntegerPow<9>(x)*IntegerPow<4>(y) - 1716*IntegerPow<7>(x)*IntegerPow<6>(y) + 1287*IntegerPow<5>(x)*IntegerPow<8>(y) - 286*IntegerPow<3>(x)*IntegerPow<10>(y) + 13*x*IntegerPow<12>(y))*(1 - 140*IntegerPow<2>(z) + 37*IntegerPow<4>(z)*(70 - 364*IntegerPow<2>(z) + 533*IntegerPow<4>(z))))/1.048576e6;
    if constexpr(I == 476) return (3*sqrt(2787362985/(2.*pi))*(IntegerPow<14>(x) - 91*IntegerPow<12>(x)*IntegerPow<2>(y) + 1001*IntegerPow<10>(x)*IntegerPow<4>(y) - 3003*IntegerPow<8>(x)*IntegerPow<6>(y) + 3003*IntegerPow<6>(x)*IntegerPow<8>(y) - 1001*IntegerPow<4>(x)*IntegerPow<10>(y) + 91*IntegerPow<2>(x)*IntegerPow<12>(y) - IntegerPow<14>(y))*z*(-35 + 37*IntegerPow<2>(z)*(35 - 273*IntegerPow<2>(z) + 533*IntegerPow<4>(z))))/262144.;
    if constexpr(I == 477) return (sqrt(19511540895/(2.*pi))*x*(IntegerPow<2>(x) - 3*IntegerPow<2>(y))*(IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + 5*IntegerPow<4>(y))*(IntegerPow<8>(x) - 92*IntegerPow<6>(x)*IntegerPow<2>(y) + 134*IntegerPow<4>(x)*IntegerPow<4>(y) - 28*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*(-5 + 37*IntegerPow<2>(z)*(15 - 195*IntegerPow<2>(z) + 533*IntegerPow<4>(z))))/524288.;
    if constexpr(I == 478) return (3*sqrt(240642337705/pi)*(IntegerPow<16>(x) - 120*IntegerPow<14>(x)*IntegerPow<2>(y) + 1820*IntegerPow<12>(x)*IntegerPow<4>(y) - 8008*IntegerPow<10>(x)*IntegerPow<6>(y) + 12870*IntegerPow<8>(x)*IntegerPow<8>(y) - 8008*IntegerPow<6>(x)*IntegerPow<10>(y) + 1820*IntegerPow<4>(x)*IntegerPow<12>(y) - 120*IntegerPow<2>(x)*IntegerPow<14>(y) + IntegerPow<16>(y))*z*(5 - 130*IntegerPow<2>(z) + 533*IntegerPow<4>(z)))/524288.;
    if constexpr(I == 479) return (15*sqrt(2533077239/(2.*pi))*(IntegerPow<17>(x) - 136*IntegerPow<15>(x)*IntegerPow<2>(y) + 2380*IntegerPow<13>(x)*IntegerPow<4>(y) - 12376*IntegerPow<11>(x)*IntegerPow<6>(y) + 24310*IntegerPow<9>(x)*IntegerPow<8>(y) - 19448*IntegerPow<7>(x)*IntegerPow<10>(y) + 6188*IntegerPow<5>(x)*IntegerPow<12>(y) - 680*IntegerPow<3>(x)*IntegerPow<14>(y) + 17*x*IntegerPow<16>(y))*(1 - 78*IntegerPow<2>(z) + 533*IntegerPow<4>(z)))/524288.;
    if constexpr(I == 480) return (5*sqrt(98790012321/(2.*pi))*(IntegerPow<18>(x) - 153*IntegerPow<16>(x)*IntegerPow<2>(y) + 3060*IntegerPow<14>(x)*IntegerPow<4>(y) - 18564*IntegerPow<12>(x)*IntegerPow<6>(y) + 43758*IntegerPow<10>(x)*IntegerPow<8>(y) - 43758*IntegerPow<8>(x)*IntegerPow<10>(y) + 18564*IntegerPow<6>(x)*IntegerPow<12>(y) - 3060*IntegerPow<4>(x)*IntegerPow<14>(y) + 153*IntegerPow<2>(x)*IntegerPow<16>(y) - IntegerPow<18>(y))*z*(-3 + 41*IntegerPow<2>(z)))/262144.;
    if constexpr(I == 481) return (3*sqrt(164650020535/pi)*(IntegerPow<19>(x) - 171*IntegerPow<17>(x)*IntegerPow<2>(y) + 3876*IntegerPow<15>(x)*IntegerPow<4>(y) - 27132*IntegerPow<13>(x)*IntegerPow<6>(y) + 75582*IntegerPow<11>(x)*IntegerPow<8>(y) - 92378*IntegerPow<9>(x)*IntegerPow<10>(y) + 50388*IntegerPow<7>(x)*IntegerPow<12>(y) - 11628*IntegerPow<5>(x)*IntegerPow<14>(y) + 969*IntegerPow<3>(x)*IntegerPow<16>(y) - 19*x*IntegerPow<18>(y))*(-1 + 41*IntegerPow<2>(z)))/1.048576e6;
    if constexpr(I == 482) return (3*sqrt(6750650841935/(2.*pi))*(IntegerPow<20>(x) - 190*IntegerPow<18>(x)*IntegerPow<2>(y) + 4845*IntegerPow<16>(x)*IntegerPow<4>(y) - 38760*IntegerPow<14>(x)*IntegerPow<6>(y) + 125970*IntegerPow<12>(x)*IntegerPow<8>(y) - 184756*IntegerPow<10>(x)*IntegerPow<10>(y) + 125970*IntegerPow<8>(x)*IntegerPow<12>(y) - 38760*IntegerPow<6>(x)*IntegerPow<14>(y) + 4845*IntegerPow<4>(x)*IntegerPow<16>(y) - 190*IntegerPow<2>(x)*IntegerPow<18>(y) + IntegerPow<20>(y))*z)/524288.;
    if constexpr(I == 483) return (sqrt(2893136075115/pi)*x*(IntegerPow<2>(x) - 3*IntegerPow<2>(y))*(IntegerPow<6>(x) - 21*IntegerPow<4>(x)*IntegerPow<2>(y) + 35*IntegerPow<2>(x)*IntegerPow<4>(y) - 7*IntegerPow<6>(y))*(IntegerPow<12>(x) - 186*IntegerPow<10>(x)*IntegerPow<2>(y) + 1423*IntegerPow<8>(x)*IntegerPow<4>(y) - 1772*IntegerPow<6>(x)*IntegerPow<6>(y) + 655*IntegerPow<4>(x)*IntegerPow<8>(y) - 58*IntegerPow<2>(x)*IntegerPow<10>(y) + IntegerPow<12>(y)))/1.048576e6;
    if constexpr(I == 484) return (15*sqrt(52602474093/pi)*x*y*(IntegerPow<10>(x) - 55*IntegerPow<8>(x)*IntegerPow<2>(y) + 330*IntegerPow<6>(x)*IntegerPow<4>(y) - 462*IntegerPow<4>(x)*IntegerPow<6>(y) + 165*IntegerPow<2>(x)*IntegerPow<8>(y) - 11*IntegerPow<10>(y))*(11*IntegerPow<10>(x) - 165*IntegerPow<8>(x)*IntegerPow<2>(y) + 462*IntegerPow<6>(x)*IntegerPow<4>(y) - 330*IntegerPow<4>(x)*IntegerPow<6>(y) + 55*IntegerPow<2>(x)*IntegerPow<8>(y) - IntegerPow<10>(y)))/1.048576e6;
    if constexpr(I == 485) return (15*sqrt(578627215023/pi)*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*(-7*IntegerPow<6>(x) + 35*IntegerPow<4>(x)*IntegerPow<2>(y) - 21*IntegerPow<2>(x)*IntegerPow<4>(y) + IntegerPow<6>(y))*(IntegerPow<12>(x) - 58*IntegerPow<10>(x)*IntegerPow<2>(y) + 655*IntegerPow<8>(x)*IntegerPow<4>(y) - 1772*IntegerPow<6>(x)*IntegerPow<6>(y) + 1423*IntegerPow<4>(x)*IntegerPow<8>(y) - 186*IntegerPow<2>(x)*IntegerPow<10>(y) + IntegerPow<12>(y))*z)/1.048576e6;
    if constexpr(I == 486) return (15*sqrt(13456446861/(2.*pi))*x*(x - y)*y*(x + y)*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(IntegerPow<4>(x) - 4*IntegerPow<3>(x)*y - 14*IntegerPow<2>(x)*IntegerPow<2>(y) - 4*x*IntegerPow<3>(y) + IntegerPow<4>(y))*(IntegerPow<4>(x) + 4*IntegerPow<3>(x)*y - 14*IntegerPow<2>(x)*IntegerPow<2>(y) + 4*x*IntegerPow<3>(y) + IntegerPow<4>(y))*(IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + 5*IntegerPow<4>(y))*(-1 + 43*IntegerPow<2>(z)))/262144.;
    if constexpr(I == 487) return (-15*sqrt(94195128027/pi)*y*(-19*IntegerPow<18>(x) + 969*IntegerPow<16>(x)*IntegerPow<2>(y) - 11628*IntegerPow<14>(x)*IntegerPow<4>(y) + 50388*IntegerPow<12>(x)*IntegerPow<6>(y) - 92378*IntegerPow<10>(x)*IntegerPow<8>(y) + 75582*IntegerPow<8>(x)*IntegerPow<10>(y) - 27132*IntegerPow<6>(x)*IntegerPow<12>(y) + 3876*IntegerPow<4>(x)*IntegerPow<14>(y) - 171*IntegerPow<2>(x)*IntegerPow<16>(y) + IntegerPow<18>(y))*z*(-3 + 43*IntegerPow<2>(z)))/1.048576e6;
    if constexpr(I == 488) return (15*sqrt(2297442147/pi)*x*y*(IntegerPow<2>(x) - 3*IntegerPow<2>(y))*(3*IntegerPow<2>(x) - IntegerPow<2>(y))*(IntegerPow<6>(x) - 33*IntegerPow<4>(x)*IntegerPow<2>(y) + 27*IntegerPow<2>(x)*IntegerPow<4>(y) - 3*IntegerPow<6>(y))*(3*IntegerPow<6>(x) - 27*IntegerPow<4>(x)*IntegerPow<2>(y) + 33*IntegerPow<2>(x)*IntegerPow<4>(y) - IntegerPow<6>(y))*(3 + 41*IntegerPow<2>(z)*(-6 + 43*IntegerPow<2>(z))))/1.048576e6;
    if constexpr(I == 489) return (15*sqrt(2297442147/(2.*pi))*y*(17*IntegerPow<16>(x) - 680*IntegerPow<14>(x)*IntegerPow<2>(y) + 6188*IntegerPow<12>(x)*IntegerPow<4>(y) - 19448*IntegerPow<10>(x)*IntegerPow<6>(y) + 24310*IntegerPow<8>(x)*IntegerPow<8>(y) - 12376*IntegerPow<6>(x)*IntegerPow<10>(y) + 2380*IntegerPow<4>(x)*IntegerPow<12>(y) - 136*IntegerPow<2>(x)*IntegerPow<14>(y) + IntegerPow<16>(y))*z*(15 - 410*IntegerPow<2>(z) + 1763*IntegerPow<4>(z)))/524288.;
    if constexpr(I == 490) return (-15*sqrt(176726319/pi)*y*(-IntegerPow<15>(x) + 35*IntegerPow<13>(x)*IntegerPow<2>(y) - 273*IntegerPow<11>(x)*IntegerPow<4>(y) + 715*IntegerPow<9>(x)*IntegerPow<6>(y) - 715*IntegerPow<7>(x)*IntegerPow<8>(y) + 273*IntegerPow<5>(x)*IntegerPow<10>(y) - 35*IntegerPow<3>(x)*IntegerPow<12>(y) + x*IntegerPow<14>(y))*(-5 + 585*IntegerPow<2>(z) - 7995*IntegerPow<4>(z) + 22919*IntegerPow<6>(z)))/65536.;
    if constexpr(I == 491) return (-15*sqrt(479685723/(2.*pi))*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(IntegerPow<8>(x) - 28*IntegerPow<6>(x)*IntegerPow<2>(y) + 134*IntegerPow<4>(x)*IntegerPow<4>(y) - 92*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*z*(-35 + 13*IntegerPow<2>(z)*(105 - 861*IntegerPow<2>(z) + 1763*IntegerPow<4>(z))))/524288.;
    if constexpr(I == 492) return (15*sqrt(12964479/pi)*x*y*(IntegerPow<6>(x) - 21*IntegerPow<4>(x)*IntegerPow<2>(y) + 35*IntegerPow<2>(x)*IntegerPow<4>(y) - 7*IntegerPow<6>(y))*(7*IntegerPow<6>(x) - 35*IntegerPow<4>(x)*IntegerPow<2>(y) + 21*IntegerPow<2>(x)*IntegerPow<4>(y) - IntegerPow<6>(y))*(35 + 37*IntegerPow<2>(z)*(-140 + 2730*IntegerPow<2>(z) + 533*IntegerPow<4>(z)*(-28 + 43*IntegerPow<2>(z)))))/1.048576e6;
    if constexpr(I == 493) return (15*sqrt(12964479/pi)*y*(13*IntegerPow<12>(x) - 286*IntegerPow<10>(x)*IntegerPow<2>(y) + 1287*IntegerPow<8>(x)*IntegerPow<4>(y) - 1716*IntegerPow<6>(x)*IntegerPow<6>(y) + 715*IntegerPow<4>(x)*IntegerPow<8>(y) - 78*IntegerPow<2>(x)*IntegerPow<10>(y) + IntegerPow<12>(y))*z*(315 + 37*IntegerPow<2>(z)*(-420 + 4914*IntegerPow<2>(z) + 533*IntegerPow<4>(z)*(-36 + 43*IntegerPow<2>(z)))))/1.048576e6;
    if constexpr(I == 494) return (15*sqrt(90751353/(2.*pi))*x*(x - y)*y*(x + y)*(IntegerPow<2>(x) - 3*IntegerPow<2>(y))*(3*IntegerPow<2>(x) - IntegerPow<2>(y))*(IntegerPow<2>(x) - 4*x*y + IntegerPow<2>(y))*(IntegerPow<2>(x) + 4*x*y + IntegerPow<2>(y))*(-9 + 1575*IntegerPow<2>(z) + 37*IntegerPow<4>(z)*(-1050 + 8190*IntegerPow<2>(z) + 533*IntegerPow<4>(z)*(-45 + 43*IntegerPow<2>(z)))))/262144.;
    if constexpr(I == 495) return (-15*sqrt(140252091/pi)*y*(-11*IntegerPow<10>(x) + 165*IntegerPow<8>(x)*IntegerPow<2>(y) - 462*IntegerPow<6>(x)*IntegerPow<4>(y) + 330*IntegerPow<4>(x)*IntegerPow<6>(y) - 55*IntegerPow<2>(x)*IntegerPow<8>(y) + IntegerPow<10>(y))*z*(-99 + 5775*IntegerPow<2>(z) - 85470*IntegerPow<4>(z) + 481*IntegerPow<6>(z)*(990 + 41*IntegerPow<2>(z)*(-55 + 43*IntegerPow<2>(z)))))/1.048576e6;
    if constexpr(I == 496) return (15*sqrt(1542773001/pi)*x*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + 5*IntegerPow<4>(y))*(3 - 594*IntegerPow<2>(z) + 17325*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(-4620 + 19305*IntegerPow<2>(z) + 533*IntegerPow<4>(z)*(-66 + 43*IntegerPow<2>(z)))))/1.048576e6;
    if constexpr(I == 497) return (15*sqrt(20056049013/(2.*pi))*y*(9*IntegerPow<8>(x) - 84*IntegerPow<6>(x)*IntegerPow<2>(y) + 126*IntegerPow<4>(x)*IntegerPow<4>(y) - 36*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*z*(3 - 198*IntegerPow<2>(z) + 3465*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(-660 + 2145*IntegerPow<2>(z) - 3198*IntegerPow<4>(z) + 1763*IntegerPow<6>(z))))/262144.;
    if constexpr(I == 498) return (-15*sqrt(92424189/pi)*y*(-IntegerPow<7>(x) + 7*IntegerPow<5>(x)*IntegerPow<2>(y) - 7*IntegerPow<3>(x)*IntegerPow<4>(y) + x*IntegerPow<6>(y))*(-3 + 31*IntegerPow<2>(z)*(21 - 693*IntegerPow<2>(z) + 8085*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(-1155 + 3003*IntegerPow<2>(z) - 3731*IntegerPow<4>(z) + 1763*IntegerPow<6>(z)))))/65536.;
    if constexpr(I == 499) return (-15*sqrt(92424189/(2.*pi))*y*(-7*IntegerPow<6>(x) + 35*IntegerPow<4>(x)*IntegerPow<2>(y) - 21*IntegerPow<2>(x)*IntegerPow<4>(y) + IntegerPow<6>(y))*z*(-45 + 31*IntegerPow<2>(z)*(105 - 2079*IntegerPow<2>(z) + 17325*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(-1925 + 4095*IntegerPow<2>(z) - 4305*IntegerPow<4>(z) + 1763*IntegerPow<6>(z)))))/262144.;
    if constexpr(I == 500) return (15*sqrt(3187041/(2.*pi))*y*(3*IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 3*x*IntegerPow<4>(y))*(45 + 29*IntegerPow<2>(z)*(-360 + 31*IntegerPow<2>(z)*(420 - 5544*IntegerPow<2>(z) + 34650*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(-3080 + 5460*IntegerPow<2>(z) - 4920*IntegerPow<4>(z) + 1763*IntegerPow<6>(z))))))/524288.;
    if constexpr(I == 501) return (15*sqrt(1312311/(2.*pi))*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*z*(765 - 59160*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(1428 - 13464*IntegerPow<2>(z) + 65450*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(-4760 + 7140*IntegerPow<2>(z) - 5576*IntegerPow<4>(z) + 1763*IntegerPow<6>(z)))))/524288.;
    if constexpr(I == 502) return (-15*sqrt(437437/pi)*x*y*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(-85 + 20655*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(-9180 + 31*IntegerPow<2>(z)*(4284 - 30294*IntegerPow<2>(z) + 117810*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(-7140 + 9180*IntegerPow<2>(z) - 6273*IntegerPow<4>(z) + 1763*IntegerPow<6>(z))))))/262144.;
    if constexpr(I == 503) return (-15*sqrt(1771/(2.*pi))*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*z*(-20995 + 1700595*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(-453492 + 31*IntegerPow<2>(z)*(151164 - 831402*IntegerPow<2>(z) + 2645370*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(-135660 + 151164*IntegerPow<2>(z) - 91143*IntegerPow<4>(z) + 22919*IntegerPow<6>(z))))))/524288.;
    if constexpr(I == 504) return (3*sqrt(8855/(2.*pi))*x*y*(4199 + 5*IntegerPow<2>(z)*(-209950 + 3*IntegerPow<2>(z)*(2834325 - 43837560*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(377910 - 1662804*IntegerPow<2>(z) + 4408950*IntegerPow<4>(z) - 7170600*IntegerPow<6>(z) + 481*IntegerPow<8>(z)*(14535 - 7790*IntegerPow<2>(z) + 1763*IntegerPow<4>(z)))))))/524288.;
    if constexpr(I == 505) return (3*sqrt(1265/pi)*y*z*(88179 + 5*IntegerPow<2>(z)*(-1469650 + 35712495*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(-4534920 + 31*IntegerPow<2>(z)*(881790 - 3174444*IntegerPow<2>(z) + 7122150*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(-271320 + 13*IntegerPow<2>(z)*(17955 - 8610*IntegerPow<2>(z) + 1763*IntegerPow<4>(z))))))))/524288.;
    if constexpr(I == 506) return (3*sqrt(5/pi)*(-88179 + 23*IntegerPow<2>(z)*(969969 + 5*IntegerPow<2>(z)*(-8083075 + 3*IntegerPow<2>(z)*(43648605 - 361659870*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(1939938 - 5819814*IntegerPow<2>(z) + 11191950*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(-373065 + 13*IntegerPow<2>(z)*(21945 - 9471*IntegerPow<2>(z) + 1763*IntegerPow<4>(z)))))))))/1.048576e6;
    if constexpr(I == 507) return (3*sqrt(1265/pi)*x*z*(88179 + 5*IntegerPow<2>(z)*(-1469650 + 35712495*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(-4534920 + 31*IntegerPow<2>(z)*(881790 - 3174444*IntegerPow<2>(z) + 7122150*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(-271320 + 13*IntegerPow<2>(z)*(17955 - 8610*IntegerPow<2>(z) + 1763*IntegerPow<4>(z))))))))/524288.;
    if constexpr(I == 508) return (-3*sqrt(8855/(2.*pi))*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(4199 + 5*IntegerPow<2>(z)*(-209950 + 3*IntegerPow<2>(z)*(2834325 - 43837560*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(377910 - 1662804*IntegerPow<2>(z) + 4408950*IntegerPow<4>(z) - 7170600*IntegerPow<6>(z) + 481*IntegerPow<8>(z)*(14535 - 7790*IntegerPow<2>(z) + 1763*IntegerPow<4>(z)))))))/1.048576e6;
    if constexpr(I == 509) return (15*sqrt(1771/(2.*pi))*(IntegerPow<3>(x) - 3*x*IntegerPow<2>(y))*z*(-20995 + 1700595*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(-453492 + 31*IntegerPow<2>(z)*(151164 - 831402*IntegerPow<2>(z) + 2645370*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(-135660 + 151164*IntegerPow<2>(z) - 91143*IntegerPow<4>(z) + 22919*IntegerPow<6>(z))))))/524288.;
    if constexpr(I == 510) return (15*sqrt(437437/pi)*(8*IntegerPow<4>(y) + 8*IntegerPow<2>(y)*(-1 + IntegerPow<2>(z)) + IntegerPow<2>(-1 + IntegerPow<2>(z)))*(-85 + 20655*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(-9180 + 31*IntegerPow<2>(z)*(4284 - 30294*IntegerPow<2>(z) + 117810*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(-7140 + 9180*IntegerPow<2>(z) - 6273*IntegerPow<4>(z) + 1763*IntegerPow<6>(z))))))/1.048576e6;
    if constexpr(I == 511) return (15*sqrt(1312311/(2.*pi))*(IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 5*x*IntegerPow<4>(y))*z*(765 - 59160*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(1428 - 13464*IntegerPow<2>(z) + 65450*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(-4760 + 7140*IntegerPow<2>(z) - 5576*IntegerPow<4>(z) + 1763*IntegerPow<6>(z)))))/524288.;
    if constexpr(I == 512) return (-15*sqrt(3187041/(2.*pi))*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(16*IntegerPow<4>(y) + 16*IntegerPow<2>(y)*(-1 + IntegerPow<2>(z)) + IntegerPow<2>(-1 + IntegerPow<2>(z)))*(45 + 29*IntegerPow<2>(z)*(-360 + 31*IntegerPow<2>(z)*(420 - 5544*IntegerPow<2>(z) + 34650*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(-3080 + 5460*IntegerPow<2>(z) - 4920*IntegerPow<4>(z) + 1763*IntegerPow<6>(z))))))/1.048576e6;
    if constexpr(I == 513) return (15*sqrt(92424189/(2.*pi))*(IntegerPow<7>(x) - 21*IntegerPow<5>(x)*IntegerPow<2>(y) + 35*IntegerPow<3>(x)*IntegerPow<4>(y) - 7*x*IntegerPow<6>(y))*z*(-45 + 31*IntegerPow<2>(z)*(105 - 2079*IntegerPow<2>(z) + 17325*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(-1925 + 4095*IntegerPow<2>(z) - 4305*IntegerPow<4>(z) + 1763*IntegerPow<6>(z)))))/262144.;
    if constexpr(I == 514) return (15*sqrt(92424189/pi)*(IntegerPow<8>(x) - 28*IntegerPow<6>(x)*IntegerPow<2>(y) + 70*IntegerPow<4>(x)*IntegerPow<4>(y) - 28*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*(-3 + 31*IntegerPow<2>(z)*(21 - 693*IntegerPow<2>(z) + 8085*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(-1155 + 3003*IntegerPow<2>(z) - 3731*IntegerPow<4>(z) + 1763*IntegerPow<6>(z)))))/524288.;
    if constexpr(I == 515) return (15*sqrt(20056049013/(2.*pi))*(IntegerPow<9>(x) - 36*IntegerPow<7>(x)*IntegerPow<2>(y) + 126*IntegerPow<5>(x)*IntegerPow<4>(y) - 84*IntegerPow<3>(x)*IntegerPow<6>(y) + 9*x*IntegerPow<8>(y))*z*(3 - 198*IntegerPow<2>(z) + 3465*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(-660 + 2145*IntegerPow<2>(z) - 3198*IntegerPow<4>(z) + 1763*IntegerPow<6>(z))))/262144.;
    if constexpr(I == 516) return (-15*sqrt(1542773001/pi)*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(256*IntegerPow<8>(y) + 512*IntegerPow<6>(y)*(-1 + IntegerPow<2>(z)) + 304*IntegerPow<4>(y)*IntegerPow<2>(-1 + IntegerPow<2>(z)) + 48*IntegerPow<2>(y)*IntegerPow<3>(-1 + IntegerPow<2>(z)) + IntegerPow<4>(-1 + IntegerPow<2>(z)))*(3 - 594*IntegerPow<2>(z) + 17325*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(-4620 + 19305*IntegerPow<2>(z) + 533*IntegerPow<4>(z)*(-66 + 43*IntegerPow<2>(z)))))/2.097152e6;
    if constexpr(I == 517) return (15*sqrt(140252091/pi)*(IntegerPow<11>(x) - 55*IntegerPow<9>(x)*IntegerPow<2>(y) + 330*IntegerPow<7>(x)*IntegerPow<4>(y) - 462*IntegerPow<5>(x)*IntegerPow<6>(y) + 165*IntegerPow<3>(x)*IntegerPow<8>(y) - 11*x*IntegerPow<10>(y))*z*(-99 + 5775*IntegerPow<2>(z) - 85470*IntegerPow<4>(z) + 481*IntegerPow<6>(z)*(990 + 41*IntegerPow<2>(z)*(-55 + 43*IntegerPow<2>(z)))))/1.048576e6;
    if constexpr(I == 518) return (15*sqrt(90751353/(2.*pi))*(8*IntegerPow<4>(y) + 8*IntegerPow<2>(y)*(-1 + IntegerPow<2>(z)) + IntegerPow<2>(-1 + IntegerPow<2>(z)))*(256*IntegerPow<8>(y) + 512*IntegerPow<6>(y)*(-1 + IntegerPow<2>(z)) + 320*IntegerPow<4>(y)*IntegerPow<2>(-1 + IntegerPow<2>(z)) + 64*IntegerPow<2>(y)*IntegerPow<3>(-1 + IntegerPow<2>(z)) + IntegerPow<4>(-1 + IntegerPow<2>(z)))*(-9 + 1575*IntegerPow<2>(z) + 37*IntegerPow<4>(z)*(-1050 + 8190*IntegerPow<2>(z) + 533*IntegerPow<4>(z)*(-45 + 43*IntegerPow<2>(z)))))/1.048576e6;
    if constexpr(I == 519) return (15*sqrt(12964479/pi)*(IntegerPow<13>(x) - 78*IntegerPow<11>(x)*IntegerPow<2>(y) + 715*IntegerPow<9>(x)*IntegerPow<4>(y) - 1716*IntegerPow<7>(x)*IntegerPow<6>(y) + 1287*IntegerPow<5>(x)*IntegerPow<8>(y) - 286*IntegerPow<3>(x)*IntegerPow<10>(y) + 13*x*IntegerPow<12>(y))*z*(315 + 37*IntegerPow<2>(z)*(-420 + 4914*IntegerPow<2>(z) + 533*IntegerPow<4>(z)*(-36 + 43*IntegerPow<2>(z)))))/1.048576e6;
    if constexpr(I == 520) return (-15*sqrt(12964479/pi)*(8192*IntegerPow<14>(y) + 28672*IntegerPow<12>(y)*(-1 + IntegerPow<2>(z)) + 39424*IntegerPow<10>(y)*IntegerPow<2>(-1 + IntegerPow<2>(z)) + 26880*IntegerPow<8>(y)*IntegerPow<3>(-1 + IntegerPow<2>(z)) + 9408*IntegerPow<6>(y)*IntegerPow<4>(-1 + IntegerPow<2>(z)) + 1568*IntegerPow<4>(y)*IntegerPow<5>(-1 + IntegerPow<2>(z)) + 98*IntegerPow<2>(y)*IntegerPow<6>(-1 + IntegerPow<2>(z)) + IntegerPow<7>(-1 + IntegerPow<2>(z)))*(35 + 37*IntegerPow<2>(z)*(-140 + 2730*IntegerPow<2>(z) + 533*IntegerPow<4>(z)*(-28 + 43*IntegerPow<2>(z)))))/2.097152e6;
    if constexpr(I == 521) return (15*sqrt(479685723/(2.*pi))*x*(IntegerPow<2>(x) - 3*IntegerPow<2>(y))*(IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + 5*IntegerPow<4>(y))*(IntegerPow<8>(x) - 92*IntegerPow<6>(x)*IntegerPow<2>(y) + 134*IntegerPow<4>(x)*IntegerPow<4>(y) - 28*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*z*(-35 + 13*IntegerPow<2>(z)*(105 - 861*IntegerPow<2>(z) + 1763*IntegerPow<4>(z))))/524288.;
    if constexpr(I == 522) return (15*sqrt(176726319/pi)*(IntegerPow<16>(x) - 120*IntegerPow<14>(x)*IntegerPow<2>(y) + 1820*IntegerPow<12>(x)*IntegerPow<4>(y) - 8008*IntegerPow<10>(x)*IntegerPow<6>(y) + 12870*IntegerPow<8>(x)*IntegerPow<8>(y) - 8008*IntegerPow<6>(x)*IntegerPow<10>(y) + 1820*IntegerPow<4>(x)*IntegerPow<12>(y) - 120*IntegerPow<2>(x)*IntegerPow<14>(y) + IntegerPow<16>(y))*(-5 + 585*IntegerPow<2>(z) - 7995*IntegerPow<4>(z) + 22919*IntegerPow<6>(z)))/1.048576e6;
    if constexpr(I == 523) return (15*sqrt(2297442147/(2.*pi))*(IntegerPow<17>(x) - 136*IntegerPow<15>(x)*IntegerPow<2>(y) + 2380*IntegerPow<13>(x)*IntegerPow<4>(y) - 12376*IntegerPow<11>(x)*IntegerPow<6>(y) + 24310*IntegerPow<9>(x)*IntegerPow<8>(y) - 19448*IntegerPow<7>(x)*IntegerPow<10>(y) + 6188*IntegerPow<5>(x)*IntegerPow<12>(y) - 680*IntegerPow<3>(x)*IntegerPow<14>(y) + 17*x*IntegerPow<16>(y))*z*(15 - 410*IntegerPow<2>(z) + 1763*IntegerPow<4>(z)))/524288.;
    if constexpr(I == 524) return (15*sqrt(2297442147/pi)*(IntegerPow<18>(x) - 153*IntegerPow<16>(x)*IntegerPow<2>(y) + 3060*IntegerPow<14>(x)*IntegerPow<4>(y) - 18564*IntegerPow<12>(x)*IntegerPow<6>(y) + 43758*IntegerPow<10>(x)*IntegerPow<8>(y) - 43758*IntegerPow<8>(x)*IntegerPow<10>(y) + 18564*IntegerPow<6>(x)*IntegerPow<12>(y) - 3060*IntegerPow<4>(x)*IntegerPow<14>(y) + 153*IntegerPow<2>(x)*IntegerPow<16>(y) - IntegerPow<18>(y))*(3 + 41*IntegerPow<2>(z)*(-6 + 43*IntegerPow<2>(z))))/2.097152e6;
    if constexpr(I == 525) return (15*sqrt(94195128027/pi)*(IntegerPow<19>(x) - 171*IntegerPow<17>(x)*IntegerPow<2>(y) + 3876*IntegerPow<15>(x)*IntegerPow<4>(y) - 27132*IntegerPow<13>(x)*IntegerPow<6>(y) + 75582*IntegerPow<11>(x)*IntegerPow<8>(y) - 92378*IntegerPow<9>(x)*IntegerPow<10>(y) + 50388*IntegerPow<7>(x)*IntegerPow<12>(y) - 11628*IntegerPow<5>(x)*IntegerPow<14>(y) + 969*IntegerPow<3>(x)*IntegerPow<16>(y) - 19*x*IntegerPow<18>(y))*z*(-3 + 43*IntegerPow<2>(z)))/1.048576e6;
    if constexpr(I == 526) return (15*sqrt(13456446861/(2.*pi))*(IntegerPow<20>(x) - 190*IntegerPow<18>(x)*IntegerPow<2>(y) + 4845*IntegerPow<16>(x)*IntegerPow<4>(y) - 38760*IntegerPow<14>(x)*IntegerPow<6>(y) + 125970*IntegerPow<12>(x)*IntegerPow<8>(y) - 184756*IntegerPow<10>(x)*IntegerPow<10>(y) + 125970*IntegerPow<8>(x)*IntegerPow<12>(y) - 38760*IntegerPow<6>(x)*IntegerPow<14>(y) + 4845*IntegerPow<4>(x)*IntegerPow<16>(y) - 190*IntegerPow<2>(x)*IntegerPow<18>(y) + IntegerPow<20>(y))*(-1 + 43*IntegerPow<2>(z)))/1.048576e6;
    if constexpr(I == 527) return (15*sqrt(578627215023/pi)*x*(IntegerPow<2>(x) - 3*IntegerPow<2>(y))*(IntegerPow<6>(x) - 21*IntegerPow<4>(x)*IntegerPow<2>(y) + 35*IntegerPow<2>(x)*IntegerPow<4>(y) - 7*IntegerPow<6>(y))*(IntegerPow<12>(x) - 186*IntegerPow<10>(x)*IntegerPow<2>(y) + 1423*IntegerPow<8>(x)*IntegerPow<4>(y) - 1772*IntegerPow<6>(x)*IntegerPow<6>(y) + 655*IntegerPow<4>(x)*IntegerPow<8>(y) - 58*IntegerPow<2>(x)*IntegerPow<10>(y) + IntegerPow<12>(y))*z)/1.048576e6;
    if constexpr(I == 528) return (15*sqrt(52602474093/pi)*(IntegerPow<22>(x) - 231*IntegerPow<20>(x)*IntegerPow<2>(y) + 7315*IntegerPow<18>(x)*IntegerPow<4>(y) - 74613*IntegerPow<16>(x)*IntegerPow<6>(y) + 319770*IntegerPow<14>(x)*IntegerPow<8>(y) - 646646*IntegerPow<12>(x)*IntegerPow<10>(y) + 646646*IntegerPow<10>(x)*IntegerPow<12>(y) - 319770*IntegerPow<8>(x)*IntegerPow<14>(y) + 74613*IntegerPow<6>(x)*IntegerPow<16>(y) - 7315*IntegerPow<4>(x)*IntegerPow<18>(y) + 231*IntegerPow<2>(x)*IntegerPow<20>(y) - IntegerPow<22>(y)))/2.097152e6;
    if constexpr(I == 529) return (-15*sqrt(107492012277/(2.*pi))*y*(-23*IntegerPow<22>(x) + 1771*IntegerPow<20>(x)*IntegerPow<2>(y) - 33649*IntegerPow<18>(x)*IntegerPow<4>(y) + 245157*IntegerPow<16>(x)*IntegerPow<6>(y) - 817190*IntegerPow<14>(x)*IntegerPow<8>(y) + 1352078*IntegerPow<12>(x)*IntegerPow<10>(y) - 1144066*IntegerPow<10>(x)*IntegerPow<12>(y) + 490314*IntegerPow<8>(x)*IntegerPow<14>(y) - 100947*IntegerPow<6>(x)*IntegerPow<16>(y) + 8855*IntegerPow<4>(x)*IntegerPow<18>(y) - 253*IntegerPow<2>(x)*IntegerPow<20>(y) + IntegerPow<22>(y)))/2.097152e6;
    if constexpr(I == 530) return (15*sqrt(2472316282371/pi)*x*y*(IntegerPow<10>(x) - 55*IntegerPow<8>(x)*IntegerPow<2>(y) + 330*IntegerPow<6>(x)*IntegerPow<4>(y) - 462*IntegerPow<4>(x)*IntegerPow<6>(y) + 165*IntegerPow<2>(x)*IntegerPow<8>(y) - 11*IntegerPow<10>(y))*(11*IntegerPow<10>(x) - 165*IntegerPow<8>(x)*IntegerPow<2>(y) + 462*IntegerPow<6>(x)*IntegerPow<4>(y) - 330*IntegerPow<4>(x)*IntegerPow<6>(y) + 55*IntegerPow<2>(x)*IntegerPow<8>(y) - IntegerPow<10>(y))*z)/1.048576e6;
    if constexpr(I == 531) return (sqrt(12361581411855/(2.*pi))*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*(-7*IntegerPow<6>(x) + 35*IntegerPow<4>(x)*IntegerPow<2>(y) - 21*IntegerPow<2>(x)*IntegerPow<4>(y) + IntegerPow<6>(y))*(IntegerPow<12>(x) - 58*IntegerPow<10>(x)*IntegerPow<2>(y) + 655*IntegerPow<8>(x)*IntegerPow<4>(y) - 1772*IntegerPow<6>(x)*IntegerPow<6>(y) + 1423*IntegerPow<4>(x)*IntegerPow<8>(y) - 186*IntegerPow<2>(x)*IntegerPow<10>(y) + IntegerPow<12>(y))*(-1 + 45*IntegerPow<2>(z)))/2.097152e6;
    if constexpr(I == 532) return (3*sqrt(45325798510135/(2.*pi))*x*(x - y)*y*(x + y)*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(IntegerPow<4>(x) - 4*IntegerPow<3>(x)*y - 14*IntegerPow<2>(x)*IntegerPow<2>(y) - 4*x*IntegerPow<3>(y) + IntegerPow<4>(y))*(IntegerPow<4>(x) + 4*IntegerPow<3>(x)*y - 14*IntegerPow<2>(x)*IntegerPow<2>(y) + 4*x*IntegerPow<3>(y) + IntegerPow<4>(y))*(IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + 5*IntegerPow<4>(y))*z*(-1 + 15*IntegerPow<2>(z)))/262144.;
    if constexpr(I == 533) return (-3*sqrt(1054088337445/(2.*pi))*y*(-19*IntegerPow<18>(x) + 969*IntegerPow<16>(x)*IntegerPow<2>(y) - 11628*IntegerPow<14>(x)*IntegerPow<4>(y) + 50388*IntegerPow<12>(x)*IntegerPow<6>(y) - 92378*IntegerPow<10>(x)*IntegerPow<8>(y) + 75582*IntegerPow<8>(x)*IntegerPow<10>(y) - 27132*IntegerPow<6>(x)*IntegerPow<12>(y) + 3876*IntegerPow<4>(x)*IntegerPow<14>(y) - 171*IntegerPow<2>(x)*IntegerPow<16>(y) + IntegerPow<18>(y))*(1 - 86*IntegerPow<2>(z) + 645*IntegerPow<4>(z)))/2.097152e6;
    if constexpr(I == 534) return (5*sqrt(4427171017269/pi)*x*y*(IntegerPow<2>(x) - 3*IntegerPow<2>(y))*(3*IntegerPow<2>(x) - IntegerPow<2>(y))*(IntegerPow<6>(x) - 33*IntegerPow<4>(x)*IntegerPow<2>(y) + 27*IntegerPow<2>(x)*IntegerPow<4>(y) - 3*IntegerPow<6>(y))*(3*IntegerPow<6>(x) - 27*IntegerPow<4>(x)*IntegerPow<2>(y) + 33*IntegerPow<2>(x)*IntegerPow<4>(y) - IntegerPow<6>(y))*z*(3 - 86*IntegerPow<2>(z) + 387*IntegerPow<4>(z)))/1.048576e6;
    if constexpr(I == 535) return (15*sqrt(35993260303/(2.*pi))*y*(17*IntegerPow<16>(x) - 680*IntegerPow<14>(x)*IntegerPow<2>(y) + 6188*IntegerPow<12>(x)*IntegerPow<4>(y) - 19448*IntegerPow<10>(x)*IntegerPow<6>(y) + 24310*IntegerPow<8>(x)*IntegerPow<8>(y) - 12376*IntegerPow<6>(x)*IntegerPow<10>(y) + 2380*IntegerPow<4>(x)*IntegerPow<12>(y) - 136*IntegerPow<2>(x)*IntegerPow<14>(y) + IntegerPow<16>(y))*(-1 + 41*IntegerPow<2>(z)*(3 - 43*IntegerPow<2>(z) + 129*IntegerPow<4>(z))))/2.097152e6;
    if constexpr(I == 536) return (-3*sqrt(25709471645/pi)*y*(-IntegerPow<15>(x) + 35*IntegerPow<13>(x)*IntegerPow<2>(y) - 273*IntegerPow<11>(x)*IntegerPow<4>(y) + 715*IntegerPow<9>(x)*IntegerPow<6>(y) - 715*IntegerPow<7>(x)*IntegerPow<8>(y) + 273*IntegerPow<5>(x)*IntegerPow<10>(y) - 35*IntegerPow<3>(x)*IntegerPow<12>(y) + x*IntegerPow<14>(y))*z*(-35 + 41*IntegerPow<2>(z)*(35 - 301*IntegerPow<2>(z) + 645*IntegerPow<4>(z))))/65536.;
    if constexpr(I == 537) return -(sqrt(5932954995/(2.*pi))*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(IntegerPow<8>(x) - 28*IntegerPow<6>(x)*IntegerPow<2>(y) + 134*IntegerPow<4>(x)*IntegerPow<4>(y) - 92*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*(35 - 5460*IntegerPow<2>(z) + 533*IntegerPow<4>(z)*(210 + 43*IntegerPow<2>(z)*(-28 + 45*IntegerPow<2>(z)))))/2.097152e6;
    if constexpr(I == 538) return (3*sqrt(112726144905/pi)*x*y*(IntegerPow<6>(x) - 21*IntegerPow<4>(x)*IntegerPow<2>(y) + 35*IntegerPow<2>(x)*IntegerPow<4>(y) - 7*IntegerPow<6>(y))*(7*IntegerPow<6>(x) - 35*IntegerPow<4>(x)*IntegerPow<2>(y) + 21*IntegerPow<2>(x)*IntegerPow<4>(y) - IntegerPow<6>(y))*z*(35 - 1820*IntegerPow<2>(z) + 533*IntegerPow<4>(z)*(42 + 43*IntegerPow<2>(z)*(-4 + 5*IntegerPow<2>(z)))))/1.048576e6;
    if constexpr(I == 539) return (15*sqrt(609330513/(2.*pi))*y*(13*IntegerPow<12>(x) - 286*IntegerPow<10>(x)*IntegerPow<2>(y) + 1287*IntegerPow<8>(x)*IntegerPow<4>(y) - 1716*IntegerPow<6>(x)*IntegerPow<6>(y) + 715*IntegerPow<4>(x)*IntegerPow<8>(y) - 78*IntegerPow<2>(x)*IntegerPow<10>(y) + IntegerPow<12>(y))*(-7 + 37*IntegerPow<2>(z)*(35 - 910*IntegerPow<2>(z) + 533*IntegerPow<4>(z)*(14 - 43*IntegerPow<2>(z) + 43*IntegerPow<4>(z)))))/2.097152e6;
    if constexpr(I == 540) return (5*sqrt(55393683/(2.*pi))*x*(x - y)*y*(x + y)*(IntegerPow<2>(x) - 3*IntegerPow<2>(y))*(3*IntegerPow<2>(x) - IntegerPow<2>(y))*(IntegerPow<2>(x) - 4*x*y + IntegerPow<2>(y))*(IntegerPow<2>(x) + 4*x*y + IntegerPow<2>(y))*z*(-693 + 37*IntegerPow<2>(z)*(1155 - 18018*IntegerPow<2>(z) + 533*IntegerPow<4>(z)*(198 - 473*IntegerPow<2>(z) + 387*IntegerPow<4>(z)))))/262144.;
    if constexpr(I == 541) return (-3*sqrt(646259635/(2.*pi))*y*(-11*IntegerPow<10>(x) + 165*IntegerPow<8>(x)*IntegerPow<2>(y) - 462*IntegerPow<6>(x)*IntegerPow<4>(y) + 330*IntegerPow<4>(x)*IntegerPow<6>(y) - 55*IntegerPow<2>(x)*IntegerPow<8>(y) + IntegerPow<10>(y))*(33 - 6930*IntegerPow<2>(z) + 37*IntegerPow<4>(z)*(5775 - 60060*IntegerPow<2>(z) + 533*IntegerPow<4>(z)*(495 - 946*IntegerPow<2>(z) + 645*IntegerPow<4>(z)))))/2.097152e6;
    if constexpr(I == 542) return (3*sqrt(142823379335/pi)*x*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + 5*IntegerPow<4>(y))*z*(33 - 2310*IntegerPow<2>(z) + 37*IntegerPow<4>(z)*(1155 - 8580*IntegerPow<2>(z) + 41*IntegerPow<4>(z)*(715 - 1118*IntegerPow<2>(z) + 645*IntegerPow<4>(z)))))/1.048576e6;
    if constexpr(I == 543) return (sqrt(673310216865/(2.*pi))*y*(9*IntegerPow<8>(x) - 84*IntegerPow<6>(x)*IntegerPow<2>(y) + 126*IntegerPow<4>(x)*IntegerPow<4>(y) - 36*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*(-3 + 693*IntegerPow<2>(z) - 24255*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(8085 - 45045*IntegerPow<2>(z) + 41*IntegerPow<4>(z)*(3003 + 43*IntegerPow<2>(z)*(-91 + 45*IntegerPow<2>(z))))))/2.097152e6;
    if constexpr(I == 544) return (-15*sqrt(44887347791/pi)*y*(-IntegerPow<7>(x) + 7*IntegerPow<5>(x)*IntegerPow<2>(y) - 7*IntegerPow<3>(x)*IntegerPow<4>(y) + x*IntegerPow<6>(y))*z*(-3 + 231*IntegerPow<2>(z) - 4851*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(1155 - 5005*IntegerPow<2>(z) + 41*IntegerPow<4>(z)*(273 + 43*IntegerPow<2>(z)*(-7 + 3*IntegerPow<2>(z))))))/65536.;
    if constexpr(I == 545) return (-15*sqrt(1447978961/pi)*y*(-7*IntegerPow<6>(x) + 35*IntegerPow<4>(x)*IntegerPow<2>(y) - 21*IntegerPow<2>(x)*IntegerPow<4>(y) + IntegerPow<6>(y))*(3 + 31*IntegerPow<2>(z)*(-24 + 924*IntegerPow<2>(z) - 12936*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(2310 - 8008*IntegerPow<2>(z) + 41*IntegerPow<4>(z)*(364 + 43*IntegerPow<2>(z)*(-8 + 3*IntegerPow<2>(z)))))))/2.097152e6;
    if constexpr(I == 546) return (sqrt(1277628495/(2.*pi))*y*(3*IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 3*x*IntegerPow<4>(y))*z*(765 + 31*IntegerPow<2>(z)*(-2040 + 47124*IntegerPow<2>(z) - 471240*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(65450 - 185640*IntegerPow<2>(z) + 41*IntegerPow<4>(z)*(7140 - 5848*IntegerPow<2>(z) + 1935*IntegerPow<4>(z))))))/524288.;
    if constexpr(I == 547) return (3*sqrt(44056155/pi)*y*(5*IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + IntegerPow<4>(y))*(-85 + 29*IntegerPow<2>(z)*(765 + 31*IntegerPow<2>(z)*(-1020 + 15708*IntegerPow<2>(z) - 117810*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(13090 - 30940*IntegerPow<2>(z) + 41820*IntegerPow<4>(z) - 29971*IntegerPow<6>(z) + 8815*IntegerPow<8>(z))))))/2.097152e6;
    if constexpr(I == 548) return (-3*sqrt(16231215/pi)*x*y*z*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(-1615 + 29*IntegerPow<2>(z)*(4845 + 31*IntegerPow<2>(z)*(-3876 + 42636*IntegerPow<2>(z) - 248710*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(22610 - 45220*IntegerPow<2>(z) + 52972*IntegerPow<4>(z) - 33497*IntegerPow<6>(z) + 8815*IntegerPow<8>(z))))))/262144.;
    if constexpr(I == 549) return (-5*sqrt(1082081/pi)*y*(-3*IntegerPow<2>(x) + IntegerPow<2>(y))*(323 - 87210*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(43605 + 31*IntegerPow<2>(z)*(-23256 + 191862*IntegerPow<2>(z) - 895356*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(67830 - 116280*IntegerPow<2>(z) + 41*IntegerPow<4>(z)*(2907 + 43*IntegerPow<2>(z)*(-38 + 9*IntegerPow<2>(z))))))))/2.097152e6;
    if constexpr(I == 550) return (5*sqrt(35673/(2.*pi))*x*y*z*(29393 - 2645370*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(793611 + 31*IntegerPow<2>(z)*(-302328 + 1939938*IntegerPow<2>(z) - 7407036*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(474810 - 705432*IntegerPow<2>(z) + 533*IntegerPow<4>(z)*(1197 - 602*IntegerPow<2>(z) + 129*IntegerPow<4>(z)))))))/524288.;
    if constexpr(I == 551) return (sqrt(3243/pi)*y*(-29393 + 8083075*IntegerPow<2>(z) + 15*IntegerPow<4>(z)*(-24249225 + 421936515*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(-4157010 + 21339318*IntegerPow<2>(z) - 67897830*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(3730650 - 4849845*IntegerPow<2>(z) + 533*IntegerPow<4>(z)*(7315 - 3311*IntegerPow<2>(z) + 645*IntegerPow<4>(z)))))))/2.097152e6;
    if constexpr(I == 552) return (sqrt(47/pi)*z*(-2028117 + 5*IntegerPow<2>(z)*(37182145 + 3*IntegerPow<2>(z)*(-334639305 + 4159088505*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(-31870410 + 133855722*IntegerPow<2>(z) - 360380790*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(17160990 - 19684665*IntegerPow<2>(z) + 533*IntegerPow<4>(z)*(26565 + 43*IntegerPow<2>(z)*(-253 + 45*IntegerPow<2>(z)))))))))/1.048576e6;
    if constexpr(I == 553) return (sqrt(3243/pi)*x*(-29393 + 8083075*IntegerPow<2>(z) + 15*IntegerPow<4>(z)*(-24249225 + 421936515*IntegerPow<2>(z) + 899*IntegerPow<4>(z)*(-4157010 + 21339318*IntegerPow<2>(z) - 67897830*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(3730650 - 4849845*IntegerPow<2>(z) + 533*IntegerPow<4>(z)*(7315 - 3311*IntegerPow<2>(z) + 645*IntegerPow<4>(z)))))))/2.097152e6;
    if constexpr(I == 554) return (-5*sqrt(35673/(2.*pi))*z*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(29393 - 2645370*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(793611 + 31*IntegerPow<2>(z)*(-302328 + 1939938*IntegerPow<2>(z) - 7407036*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(474810 - 705432*IntegerPow<2>(z) + 533*IntegerPow<4>(z)*(1197 - 602*IntegerPow<2>(z) + 129*IntegerPow<4>(z)))))))/1.048576e6;
    if constexpr(I == 555) return (5*sqrt(1082081/pi)*(IntegerPow<3>(x) - 3*x*IntegerPow<2>(y))*(323 - 87210*IntegerPow<2>(z) + 87*IntegerPow<4>(z)*(43605 + 31*IntegerPow<2>(z)*(-23256 + 191862*IntegerPow<2>(z) - 895356*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(67830 - 116280*IntegerPow<2>(z) + 41*IntegerPow<4>(z)*(2907 + 43*IntegerPow<2>(z)*(-38 + 9*IntegerPow<2>(z))))))))/2.097152e6;
    if constexpr(I == 556) return (3*sqrt(16231215/pi)*z*(8*IntegerPow<4>(y) + 8*IntegerPow<2>(y)*(-1 + IntegerPow<2>(z)) + IntegerPow<2>(-1 + IntegerPow<2>(z)))*(-1615 + 29*IntegerPow<2>(z)*(4845 + 31*IntegerPow<2>(z)*(-3876 + 42636*IntegerPow<2>(z) - 248710*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(22610 - 45220*IntegerPow<2>(z) + 52972*IntegerPow<4>(z) - 33497*IntegerPow<6>(z) + 8815*IntegerPow<8>(z))))))/1.048576e6;
    if constexpr(I == 557) return (3*sqrt(44056155/pi)*(IntegerPow<5>(x) - 10*IntegerPow<3>(x)*IntegerPow<2>(y) + 5*x*IntegerPow<4>(y))*(-85 + 29*IntegerPow<2>(z)*(765 + 31*IntegerPow<2>(z)*(-1020 + 15708*IntegerPow<2>(z) - 117810*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(13090 - 30940*IntegerPow<2>(z) + 41820*IntegerPow<4>(z) - 29971*IntegerPow<6>(z) + 8815*IntegerPow<8>(z))))))/2.097152e6;
    if constexpr(I == 558) return -(sqrt(1277628495/(2.*pi))*z*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(16*IntegerPow<4>(y) + 16*IntegerPow<2>(y)*(-1 + IntegerPow<2>(z)) + IntegerPow<2>(-1 + IntegerPow<2>(z)))*(765 + 31*IntegerPow<2>(z)*(-2040 + 47124*IntegerPow<2>(z) - 471240*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(65450 - 185640*IntegerPow<2>(z) + 41*IntegerPow<4>(z)*(7140 - 5848*IntegerPow<2>(z) + 1935*IntegerPow<4>(z))))))/1.048576e6;
    if constexpr(I == 559) return (15*sqrt(1447978961/pi)*(IntegerPow<7>(x) - 21*IntegerPow<5>(x)*IntegerPow<2>(y) + 35*IntegerPow<3>(x)*IntegerPow<4>(y) - 7*x*IntegerPow<6>(y))*(3 + 31*IntegerPow<2>(z)*(-24 + 924*IntegerPow<2>(z) - 12936*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(2310 - 8008*IntegerPow<2>(z) + 41*IntegerPow<4>(z)*(364 + 43*IntegerPow<2>(z)*(-8 + 3*IntegerPow<2>(z)))))))/2.097152e6;
    if constexpr(I == 560) return (15*sqrt(44887347791/pi)*(IntegerPow<8>(x) - 28*IntegerPow<6>(x)*IntegerPow<2>(y) + 70*IntegerPow<4>(x)*IntegerPow<4>(y) - 28*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*z*(-3 + 231*IntegerPow<2>(z) - 4851*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(1155 - 5005*IntegerPow<2>(z) + 41*IntegerPow<4>(z)*(273 + 43*IntegerPow<2>(z)*(-7 + 3*IntegerPow<2>(z))))))/524288.;
    if constexpr(I == 561) return (sqrt(673310216865/(2.*pi))*(IntegerPow<9>(x) - 36*IntegerPow<7>(x)*IntegerPow<2>(y) + 126*IntegerPow<5>(x)*IntegerPow<4>(y) - 84*IntegerPow<3>(x)*IntegerPow<6>(y) + 9*x*IntegerPow<8>(y))*(-3 + 693*IntegerPow<2>(z) - 24255*IntegerPow<4>(z) + 37*IntegerPow<6>(z)*(8085 - 45045*IntegerPow<2>(z) + 41*IntegerPow<4>(z)*(3003 + 43*IntegerPow<2>(z)*(-91 + 45*IntegerPow<2>(z))))))/2.097152e6;
    if constexpr(I == 562) return (-3*sqrt(142823379335/pi)*z*(-1 + 2*IntegerPow<2>(y) + IntegerPow<2>(z))*(256*IntegerPow<8>(y) + 512*IntegerPow<6>(y)*(-1 + IntegerPow<2>(z)) + 304*IntegerPow<4>(y)*IntegerPow<2>(-1 + IntegerPow<2>(z)) + 48*IntegerPow<2>(y)*IntegerPow<3>(-1 + IntegerPow<2>(z)) + IntegerPow<4>(-1 + IntegerPow<2>(z)))*(33 - 2310*IntegerPow<2>(z) + 37*IntegerPow<4>(z)*(1155 - 8580*IntegerPow<2>(z) + 41*IntegerPow<4>(z)*(715 - 1118*IntegerPow<2>(z) + 645*IntegerPow<4>(z)))))/2.097152e6;
    if constexpr(I == 563) return (3*sqrt(646259635/(2.*pi))*(IntegerPow<11>(x) - 55*IntegerPow<9>(x)*IntegerPow<2>(y) + 330*IntegerPow<7>(x)*IntegerPow<4>(y) - 462*IntegerPow<5>(x)*IntegerPow<6>(y) + 165*IntegerPow<3>(x)*IntegerPow<8>(y) - 11*x*IntegerPow<10>(y))*(33 - 6930*IntegerPow<2>(z) + 37*IntegerPow<4>(z)*(5775 - 60060*IntegerPow<2>(z) + 533*IntegerPow<4>(z)*(495 - 946*IntegerPow<2>(z) + 645*IntegerPow<4>(z)))))/2.097152e6;
    if constexpr(I == 564) return (5*sqrt(55393683/(2.*pi))*z*(8*IntegerPow<4>(y) + 8*IntegerPow<2>(y)*(-1 + IntegerPow<2>(z)) + IntegerPow<2>(-1 + IntegerPow<2>(z)))*(256*IntegerPow<8>(y) + 512*IntegerPow<6>(y)*(-1 + IntegerPow<2>(z)) + 320*IntegerPow<4>(y)*IntegerPow<2>(-1 + IntegerPow<2>(z)) + 64*IntegerPow<2>(y)*IntegerPow<3>(-1 + IntegerPow<2>(z)) + IntegerPow<4>(-1 + IntegerPow<2>(z)))*(-693 + 37*IntegerPow<2>(z)*(1155 - 18018*IntegerPow<2>(z) + 533*IntegerPow<4>(z)*(198 - 473*IntegerPow<2>(z) + 387*IntegerPow<4>(z)))))/1.048576e6;
    if constexpr(I == 565) return (15*sqrt(609330513/(2.*pi))*(IntegerPow<13>(x) - 78*IntegerPow<11>(x)*IntegerPow<2>(y) + 715*IntegerPow<9>(x)*IntegerPow<4>(y) - 1716*IntegerPow<7>(x)*IntegerPow<6>(y) + 1287*IntegerPow<5>(x)*IntegerPow<8>(y) - 286*IntegerPow<3>(x)*IntegerPow<10>(y) + 13*x*IntegerPow<12>(y))*(-7 + 37*IntegerPow<2>(z)*(35 - 910*IntegerPow<2>(z) + 533*IntegerPow<4>(z)*(14 - 43*IntegerPow<2>(z) + 43*IntegerPow<4>(z)))))/2.097152e6;
    if constexpr(I == 566) return (-3*sqrt(112726144905/pi)*z*(8192*IntegerPow<14>(y) + 28672*IntegerPow<12>(y)*(-1 + IntegerPow<2>(z)) + 39424*IntegerPow<10>(y)*IntegerPow<2>(-1 + IntegerPow<2>(z)) + 26880*IntegerPow<8>(y)*IntegerPow<3>(-1 + IntegerPow<2>(z)) + 9408*IntegerPow<6>(y)*IntegerPow<4>(-1 + IntegerPow<2>(z)) + 1568*IntegerPow<4>(y)*IntegerPow<5>(-1 + IntegerPow<2>(z)) + 98*IntegerPow<2>(y)*IntegerPow<6>(-1 + IntegerPow<2>(z)) + IntegerPow<7>(-1 + IntegerPow<2>(z)))*(35 - 1820*IntegerPow<2>(z) + 533*IntegerPow<4>(z)*(42 + 43*IntegerPow<2>(z)*(-4 + 5*IntegerPow<2>(z)))))/2.097152e6;
    if constexpr(I == 567) return (sqrt(5932954995/(2.*pi))*x*(IntegerPow<2>(x) - 3*IntegerPow<2>(y))*(IntegerPow<4>(x) - 10*IntegerPow<2>(x)*IntegerPow<2>(y) + 5*IntegerPow<4>(y))*(IntegerPow<8>(x) - 92*IntegerPow<6>(x)*IntegerPow<2>(y) + 134*IntegerPow<4>(x)*IntegerPow<4>(y) - 28*IntegerPow<2>(x)*IntegerPow<6>(y) + IntegerPow<8>(y))*(35 - 5460*IntegerPow<2>(z) + 533*IntegerPow<4>(z)*(210 + 43*IntegerPow<2>(z)*(-28 + 45*IntegerPow<2>(z)))))/2.097152e6;
    if constexpr(I == 568) return (3*sqrt(25709471645/pi)*(IntegerPow<16>(x) - 120*IntegerPow<14>(x)*IntegerPow<2>(y) + 1820*IntegerPow<12>(x)*IntegerPow<4>(y) - 8008*IntegerPow<10>(x)*IntegerPow<6>(y) + 12870*IntegerPow<8>(x)*IntegerPow<8>(y) - 8008*IntegerPow<6>(x)*IntegerPow<10>(y) + 1820*IntegerPow<4>(x)*IntegerPow<12>(y) - 120*IntegerPow<2>(x)*IntegerPow<14>(y) + IntegerPow<16>(y))*z*(-35 + 41*IntegerPow<2>(z)*(35 - 301*IntegerPow<2>(z) + 645*IntegerPow<4>(z))))/1.048576e6;
    if constexpr(I == 569) return (15*sqrt(35993260303/(2.*pi))*(IntegerPow<17>(x) - 136*IntegerPow<15>(x)*IntegerPow<2>(y) + 2380*IntegerPow<13>(x)*IntegerPow<4>(y) - 12376*IntegerPow<11>(x)*IntegerPow<6>(y) + 24310*IntegerPow<9>(x)*IntegerPow<8>(y) - 19448*IntegerPow<7>(x)*IntegerPow<10>(y) + 6188*IntegerPow<5>(x)*IntegerPow<12>(y) - 680*IntegerPow<3>(x)*IntegerPow<14>(y) + 17*x*IntegerPow<16>(y))*(-1 + 41*IntegerPow<2>(z)*(3 - 43*IntegerPow<2>(z) + 129*IntegerPow<4>(z))))/2.097152e6;
    if constexpr(I == 570) return (5*sqrt(4427171017269/pi)*(IntegerPow<18>(x) - 153*IntegerPow<16>(x)*IntegerPow<2>(y) + 3060*IntegerPow<14>(x)*IntegerPow<4>(y) - 18564*IntegerPow<12>(x)*IntegerPow<6>(y) + 43758*IntegerPow<10>(x)*IntegerPow<8>(y) - 43758*IntegerPow<8>(x)*IntegerPow<10>(y) + 18564*IntegerPow<6>(x)*IntegerPow<12>(y) - 3060*IntegerPow<4>(x)*IntegerPow<14>(y) + 153*IntegerPow<2>(x)*IntegerPow<16>(y) - IntegerPow<18>(y))*z*(3 - 86*IntegerPow<2>(z) + 387*IntegerPow<4>(z)))/2.097152e6;
    if constexpr(I == 571) return (3*sqrt(1054088337445/(2.*pi))*(IntegerPow<19>(x) - 171*IntegerPow<17>(x)*IntegerPow<2>(y) + 3876*IntegerPow<15>(x)*IntegerPow<4>(y) - 27132*IntegerPow<13>(x)*IntegerPow<6>(y) + 75582*IntegerPow<11>(x)*IntegerPow<8>(y) - 92378*IntegerPow<9>(x)*IntegerPow<10>(y) + 50388*IntegerPow<7>(x)*IntegerPow<12>(y) - 11628*IntegerPow<5>(x)*IntegerPow<14>(y) + 969*IntegerPow<3>(x)*IntegerPow<16>(y) - 19*x*IntegerPow<18>(y))*(1 - 86*IntegerPow<2>(z) + 645*IntegerPow<4>(z)))/2.097152e6;
    if constexpr(I == 572) return (3*sqrt(45325798510135/(2.*pi))*(IntegerPow<20>(x) - 190*IntegerPow<18>(x)*IntegerPow<2>(y) + 4845*IntegerPow<16>(x)*IntegerPow<4>(y) - 38760*IntegerPow<14>(x)*IntegerPow<6>(y) + 125970*IntegerPow<12>(x)*IntegerPow<8>(y) - 184756*IntegerPow<10>(x)*IntegerPow<10>(y) + 125970*IntegerPow<8>(x)*IntegerPow<12>(y) - 38760*IntegerPow<6>(x)*IntegerPow<14>(y) + 4845*IntegerPow<4>(x)*IntegerPow<16>(y) - 190*IntegerPow<2>(x)*IntegerPow<18>(y) + IntegerPow<20>(y))*z*(-1 + 15*IntegerPow<2>(z)))/1.048576e6;
    if constexpr(I == 573) return (sqrt(12361581411855/(2.*pi))*x*(IntegerPow<2>(x) - 3*IntegerPow<2>(y))*(IntegerPow<6>(x) - 21*IntegerPow<4>(x)*IntegerPow<2>(y) + 35*IntegerPow<2>(x)*IntegerPow<4>(y) - 7*IntegerPow<6>(y))*(IntegerPow<12>(x) - 186*IntegerPow<10>(x)*IntegerPow<2>(y) + 1423*IntegerPow<8>(x)*IntegerPow<4>(y) - 1772*IntegerPow<6>(x)*IntegerPow<6>(y) + 655*IntegerPow<4>(x)*IntegerPow<8>(y) - 58*IntegerPow<2>(x)*IntegerPow<10>(y) + IntegerPow<12>(y))*(-1 + 45*IntegerPow<2>(z)))/2.097152e6;
    if constexpr(I == 574) return (15*sqrt(2472316282371/pi)*(IntegerPow<22>(x) - 231*IntegerPow<20>(x)*IntegerPow<2>(y) + 7315*IntegerPow<18>(x)*IntegerPow<4>(y) - 74613*IntegerPow<16>(x)*IntegerPow<6>(y) + 319770*IntegerPow<14>(x)*IntegerPow<8>(y) - 646646*IntegerPow<12>(x)*IntegerPow<10>(y) + 646646*IntegerPow<10>(x)*IntegerPow<12>(y) - 319770*IntegerPow<8>(x)*IntegerPow<14>(y) + 74613*IntegerPow<6>(x)*IntegerPow<16>(y) - 7315*IntegerPow<4>(x)*IntegerPow<18>(y) + 231*IntegerPow<2>(x)*IntegerPow<20>(y) - IntegerPow<22>(y))*z)/2.097152e6;
    if constexpr(I == 575) return (15*sqrt(107492012277/(2.*pi))*(IntegerPow<23>(x) - 253*IntegerPow<21>(x)*IntegerPow<2>(y) + 8855*IntegerPow<19>(x)*IntegerPow<4>(y) - 100947*IntegerPow<17>(x)*IntegerPow<6>(y) + 490314*IntegerPow<15>(x)*IntegerPow<8>(y) - 1144066*IntegerPow<13>(x)*IntegerPow<10>(y) + 1352078*IntegerPow<11>(x)*IntegerPow<12>(y) - 817190*IntegerPow<9>(x)*IntegerPow<14>(y) + 245157*IntegerPow<7>(x)*IntegerPow<16>(y) - 33649*IntegerPow<5>(x)*IntegerPow<18>(y) + 1771*IntegerPow<3>(x)*IntegerPow<20>(y) - 23*x*IntegerPow<22>(y)))/2.097152e6;
    if constexpr(I > 575)
    {
        ExitOnError("Spherical Harmonic index i to big.");
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
  &SphericalHarmonicsXyz::Y<250>, &SphericalHarmonicsXyz::Y<251>, &SphericalHarmonicsXyz::Y<252>, &SphericalHarmonicsXyz::Y<253>, &SphericalHarmonicsXyz::Y<254>, &SphericalHarmonicsXyz::Y<255>, &SphericalHarmonicsXyz::Y<256>, &SphericalHarmonicsXyz::Y<257>, &SphericalHarmonicsXyz::Y<258>, &SphericalHarmonicsXyz::Y<259>,
  &SphericalHarmonicsXyz::Y<260>, &SphericalHarmonicsXyz::Y<261>, &SphericalHarmonicsXyz::Y<262>, &SphericalHarmonicsXyz::Y<263>, &SphericalHarmonicsXyz::Y<264>, &SphericalHarmonicsXyz::Y<265>, &SphericalHarmonicsXyz::Y<266>, &SphericalHarmonicsXyz::Y<267>, &SphericalHarmonicsXyz::Y<268>, &SphericalHarmonicsXyz::Y<269>,
  &SphericalHarmonicsXyz::Y<270>, &SphericalHarmonicsXyz::Y<271>, &SphericalHarmonicsXyz::Y<272>, &SphericalHarmonicsXyz::Y<273>, &SphericalHarmonicsXyz::Y<274>, &SphericalHarmonicsXyz::Y<275>, &SphericalHarmonicsXyz::Y<276>, &SphericalHarmonicsXyz::Y<277>, &SphericalHarmonicsXyz::Y<278>, &SphericalHarmonicsXyz::Y<279>,
  &SphericalHarmonicsXyz::Y<280>, &SphericalHarmonicsXyz::Y<281>, &SphericalHarmonicsXyz::Y<282>, &SphericalHarmonicsXyz::Y<283>, &SphericalHarmonicsXyz::Y<284>, &SphericalHarmonicsXyz::Y<285>, &SphericalHarmonicsXyz::Y<286>, &SphericalHarmonicsXyz::Y<287>, &SphericalHarmonicsXyz::Y<288>, &SphericalHarmonicsXyz::Y<289>,
  &SphericalHarmonicsXyz::Y<290>, &SphericalHarmonicsXyz::Y<291>, &SphericalHarmonicsXyz::Y<292>, &SphericalHarmonicsXyz::Y<293>, &SphericalHarmonicsXyz::Y<294>, &SphericalHarmonicsXyz::Y<295>, &SphericalHarmonicsXyz::Y<296>, &SphericalHarmonicsXyz::Y<297>, &SphericalHarmonicsXyz::Y<298>, &SphericalHarmonicsXyz::Y<299>,
  &SphericalHarmonicsXyz::Y<300>, &SphericalHarmonicsXyz::Y<301>, &SphericalHarmonicsXyz::Y<302>, &SphericalHarmonicsXyz::Y<303>, &SphericalHarmonicsXyz::Y<304>, &SphericalHarmonicsXyz::Y<305>, &SphericalHarmonicsXyz::Y<306>, &SphericalHarmonicsXyz::Y<307>, &SphericalHarmonicsXyz::Y<308>, &SphericalHarmonicsXyz::Y<309>,
  &SphericalHarmonicsXyz::Y<310>, &SphericalHarmonicsXyz::Y<311>, &SphericalHarmonicsXyz::Y<312>, &SphericalHarmonicsXyz::Y<313>, &SphericalHarmonicsXyz::Y<314>, &SphericalHarmonicsXyz::Y<315>, &SphericalHarmonicsXyz::Y<316>, &SphericalHarmonicsXyz::Y<317>, &SphericalHarmonicsXyz::Y<318>, &SphericalHarmonicsXyz::Y<319>,
  &SphericalHarmonicsXyz::Y<320>, &SphericalHarmonicsXyz::Y<321>, &SphericalHarmonicsXyz::Y<322>, &SphericalHarmonicsXyz::Y<323>, &SphericalHarmonicsXyz::Y<324>, &SphericalHarmonicsXyz::Y<325>, &SphericalHarmonicsXyz::Y<326>, &SphericalHarmonicsXyz::Y<327>, &SphericalHarmonicsXyz::Y<328>, &SphericalHarmonicsXyz::Y<329>,
  &SphericalHarmonicsXyz::Y<330>, &SphericalHarmonicsXyz::Y<331>, &SphericalHarmonicsXyz::Y<332>, &SphericalHarmonicsXyz::Y<333>, &SphericalHarmonicsXyz::Y<334>, &SphericalHarmonicsXyz::Y<335>, &SphericalHarmonicsXyz::Y<336>, &SphericalHarmonicsXyz::Y<337>, &SphericalHarmonicsXyz::Y<338>, &SphericalHarmonicsXyz::Y<339>,
  &SphericalHarmonicsXyz::Y<340>, &SphericalHarmonicsXyz::Y<341>, &SphericalHarmonicsXyz::Y<342>, &SphericalHarmonicsXyz::Y<343>, &SphericalHarmonicsXyz::Y<344>, &SphericalHarmonicsXyz::Y<345>, &SphericalHarmonicsXyz::Y<346>, &SphericalHarmonicsXyz::Y<347>, &SphericalHarmonicsXyz::Y<348>, &SphericalHarmonicsXyz::Y<349>,
  &SphericalHarmonicsXyz::Y<350>, &SphericalHarmonicsXyz::Y<351>, &SphericalHarmonicsXyz::Y<352>, &SphericalHarmonicsXyz::Y<353>, &SphericalHarmonicsXyz::Y<354>, &SphericalHarmonicsXyz::Y<355>, &SphericalHarmonicsXyz::Y<356>, &SphericalHarmonicsXyz::Y<357>, &SphericalHarmonicsXyz::Y<358>, &SphericalHarmonicsXyz::Y<359>,
  &SphericalHarmonicsXyz::Y<360>, &SphericalHarmonicsXyz::Y<361>, &SphericalHarmonicsXyz::Y<362>, &SphericalHarmonicsXyz::Y<363>, &SphericalHarmonicsXyz::Y<364>, &SphericalHarmonicsXyz::Y<365>, &SphericalHarmonicsXyz::Y<366>, &SphericalHarmonicsXyz::Y<367>, &SphericalHarmonicsXyz::Y<368>, &SphericalHarmonicsXyz::Y<369>,
  &SphericalHarmonicsXyz::Y<370>, &SphericalHarmonicsXyz::Y<371>, &SphericalHarmonicsXyz::Y<372>, &SphericalHarmonicsXyz::Y<373>, &SphericalHarmonicsXyz::Y<374>, &SphericalHarmonicsXyz::Y<375>, &SphericalHarmonicsXyz::Y<376>, &SphericalHarmonicsXyz::Y<377>, &SphericalHarmonicsXyz::Y<378>, &SphericalHarmonicsXyz::Y<379>,
  &SphericalHarmonicsXyz::Y<380>, &SphericalHarmonicsXyz::Y<381>, &SphericalHarmonicsXyz::Y<382>, &SphericalHarmonicsXyz::Y<383>, &SphericalHarmonicsXyz::Y<384>, &SphericalHarmonicsXyz::Y<385>, &SphericalHarmonicsXyz::Y<386>, &SphericalHarmonicsXyz::Y<387>, &SphericalHarmonicsXyz::Y<388>, &SphericalHarmonicsXyz::Y<389>,
  &SphericalHarmonicsXyz::Y<390>, &SphericalHarmonicsXyz::Y<391>, &SphericalHarmonicsXyz::Y<392>, &SphericalHarmonicsXyz::Y<393>, &SphericalHarmonicsXyz::Y<394>, &SphericalHarmonicsXyz::Y<395>, &SphericalHarmonicsXyz::Y<396>, &SphericalHarmonicsXyz::Y<397>, &SphericalHarmonicsXyz::Y<398>, &SphericalHarmonicsXyz::Y<399>,
  &SphericalHarmonicsXyz::Y<400>, &SphericalHarmonicsXyz::Y<401>, &SphericalHarmonicsXyz::Y<402>, &SphericalHarmonicsXyz::Y<403>, &SphericalHarmonicsXyz::Y<404>, &SphericalHarmonicsXyz::Y<405>, &SphericalHarmonicsXyz::Y<406>, &SphericalHarmonicsXyz::Y<407>, &SphericalHarmonicsXyz::Y<408>, &SphericalHarmonicsXyz::Y<409>,
  &SphericalHarmonicsXyz::Y<410>, &SphericalHarmonicsXyz::Y<411>, &SphericalHarmonicsXyz::Y<412>, &SphericalHarmonicsXyz::Y<413>, &SphericalHarmonicsXyz::Y<414>, &SphericalHarmonicsXyz::Y<415>, &SphericalHarmonicsXyz::Y<416>, &SphericalHarmonicsXyz::Y<417>, &SphericalHarmonicsXyz::Y<418>, &SphericalHarmonicsXyz::Y<419>,
  &SphericalHarmonicsXyz::Y<420>, &SphericalHarmonicsXyz::Y<421>, &SphericalHarmonicsXyz::Y<422>, &SphericalHarmonicsXyz::Y<423>, &SphericalHarmonicsXyz::Y<424>, &SphericalHarmonicsXyz::Y<425>, &SphericalHarmonicsXyz::Y<426>, &SphericalHarmonicsXyz::Y<427>, &SphericalHarmonicsXyz::Y<428>, &SphericalHarmonicsXyz::Y<429>,
  &SphericalHarmonicsXyz::Y<430>, &SphericalHarmonicsXyz::Y<431>, &SphericalHarmonicsXyz::Y<432>, &SphericalHarmonicsXyz::Y<433>, &SphericalHarmonicsXyz::Y<434>, &SphericalHarmonicsXyz::Y<435>, &SphericalHarmonicsXyz::Y<436>, &SphericalHarmonicsXyz::Y<437>, &SphericalHarmonicsXyz::Y<438>, &SphericalHarmonicsXyz::Y<439>,
  &SphericalHarmonicsXyz::Y<440>, &SphericalHarmonicsXyz::Y<441>, &SphericalHarmonicsXyz::Y<442>, &SphericalHarmonicsXyz::Y<443>, &SphericalHarmonicsXyz::Y<444>, &SphericalHarmonicsXyz::Y<445>, &SphericalHarmonicsXyz::Y<446>, &SphericalHarmonicsXyz::Y<447>, &SphericalHarmonicsXyz::Y<448>, &SphericalHarmonicsXyz::Y<449>,
  &SphericalHarmonicsXyz::Y<450>, &SphericalHarmonicsXyz::Y<451>, &SphericalHarmonicsXyz::Y<452>, &SphericalHarmonicsXyz::Y<453>, &SphericalHarmonicsXyz::Y<454>, &SphericalHarmonicsXyz::Y<455>, &SphericalHarmonicsXyz::Y<456>, &SphericalHarmonicsXyz::Y<457>, &SphericalHarmonicsXyz::Y<458>, &SphericalHarmonicsXyz::Y<459>,
  &SphericalHarmonicsXyz::Y<460>, &SphericalHarmonicsXyz::Y<461>, &SphericalHarmonicsXyz::Y<462>, &SphericalHarmonicsXyz::Y<463>, &SphericalHarmonicsXyz::Y<464>, &SphericalHarmonicsXyz::Y<465>, &SphericalHarmonicsXyz::Y<466>, &SphericalHarmonicsXyz::Y<467>, &SphericalHarmonicsXyz::Y<468>, &SphericalHarmonicsXyz::Y<469>,
  &SphericalHarmonicsXyz::Y<470>, &SphericalHarmonicsXyz::Y<471>, &SphericalHarmonicsXyz::Y<472>, &SphericalHarmonicsXyz::Y<473>, &SphericalHarmonicsXyz::Y<474>, &SphericalHarmonicsXyz::Y<475>, &SphericalHarmonicsXyz::Y<476>, &SphericalHarmonicsXyz::Y<477>, &SphericalHarmonicsXyz::Y<478>, &SphericalHarmonicsXyz::Y<479>,
  &SphericalHarmonicsXyz::Y<480>, &SphericalHarmonicsXyz::Y<481>, &SphericalHarmonicsXyz::Y<482>, &SphericalHarmonicsXyz::Y<483>, &SphericalHarmonicsXyz::Y<484>, &SphericalHarmonicsXyz::Y<485>, &SphericalHarmonicsXyz::Y<486>, &SphericalHarmonicsXyz::Y<487>, &SphericalHarmonicsXyz::Y<488>, &SphericalHarmonicsXyz::Y<489>,
  &SphericalHarmonicsXyz::Y<490>, &SphericalHarmonicsXyz::Y<491>, &SphericalHarmonicsXyz::Y<492>, &SphericalHarmonicsXyz::Y<493>, &SphericalHarmonicsXyz::Y<494>, &SphericalHarmonicsXyz::Y<495>, &SphericalHarmonicsXyz::Y<496>, &SphericalHarmonicsXyz::Y<497>, &SphericalHarmonicsXyz::Y<498>, &SphericalHarmonicsXyz::Y<499>,
  &SphericalHarmonicsXyz::Y<500>, &SphericalHarmonicsXyz::Y<501>, &SphericalHarmonicsXyz::Y<502>, &SphericalHarmonicsXyz::Y<503>, &SphericalHarmonicsXyz::Y<504>, &SphericalHarmonicsXyz::Y<505>, &SphericalHarmonicsXyz::Y<506>, &SphericalHarmonicsXyz::Y<507>, &SphericalHarmonicsXyz::Y<508>, &SphericalHarmonicsXyz::Y<509>,
  &SphericalHarmonicsXyz::Y<510>, &SphericalHarmonicsXyz::Y<511>, &SphericalHarmonicsXyz::Y<512>, &SphericalHarmonicsXyz::Y<513>, &SphericalHarmonicsXyz::Y<514>, &SphericalHarmonicsXyz::Y<515>, &SphericalHarmonicsXyz::Y<516>, &SphericalHarmonicsXyz::Y<517>, &SphericalHarmonicsXyz::Y<518>, &SphericalHarmonicsXyz::Y<519>,
  &SphericalHarmonicsXyz::Y<520>, &SphericalHarmonicsXyz::Y<521>, &SphericalHarmonicsXyz::Y<522>, &SphericalHarmonicsXyz::Y<523>, &SphericalHarmonicsXyz::Y<524>, &SphericalHarmonicsXyz::Y<525>, &SphericalHarmonicsXyz::Y<526>, &SphericalHarmonicsXyz::Y<527>, &SphericalHarmonicsXyz::Y<528>, &SphericalHarmonicsXyz::Y<529>,
  &SphericalHarmonicsXyz::Y<530>, &SphericalHarmonicsXyz::Y<531>, &SphericalHarmonicsXyz::Y<532>, &SphericalHarmonicsXyz::Y<533>, &SphericalHarmonicsXyz::Y<534>, &SphericalHarmonicsXyz::Y<535>, &SphericalHarmonicsXyz::Y<536>, &SphericalHarmonicsXyz::Y<537>, &SphericalHarmonicsXyz::Y<538>, &SphericalHarmonicsXyz::Y<539>,
  &SphericalHarmonicsXyz::Y<540>, &SphericalHarmonicsXyz::Y<541>, &SphericalHarmonicsXyz::Y<542>, &SphericalHarmonicsXyz::Y<543>, &SphericalHarmonicsXyz::Y<544>, &SphericalHarmonicsXyz::Y<545>, &SphericalHarmonicsXyz::Y<546>, &SphericalHarmonicsXyz::Y<547>, &SphericalHarmonicsXyz::Y<548>, &SphericalHarmonicsXyz::Y<549>,
  &SphericalHarmonicsXyz::Y<550>, &SphericalHarmonicsXyz::Y<551>, &SphericalHarmonicsXyz::Y<552>, &SphericalHarmonicsXyz::Y<553>, &SphericalHarmonicsXyz::Y<554>, &SphericalHarmonicsXyz::Y<555>, &SphericalHarmonicsXyz::Y<556>, &SphericalHarmonicsXyz::Y<557>, &SphericalHarmonicsXyz::Y<558>, &SphericalHarmonicsXyz::Y<559>,
  &SphericalHarmonicsXyz::Y<560>, &SphericalHarmonicsXyz::Y<561>, &SphericalHarmonicsXyz::Y<562>, &SphericalHarmonicsXyz::Y<563>, &SphericalHarmonicsXyz::Y<564>, &SphericalHarmonicsXyz::Y<565>, &SphericalHarmonicsXyz::Y<566>, &SphericalHarmonicsXyz::Y<567>, &SphericalHarmonicsXyz::Y<568>, &SphericalHarmonicsXyz::Y<569>,
  &SphericalHarmonicsXyz::Y<570>, &SphericalHarmonicsXyz::Y<571>, &SphericalHarmonicsXyz::Y<572>, &SphericalHarmonicsXyz::Y<573>, &SphericalHarmonicsXyz::Y<574>, &SphericalHarmonicsXyz::Y<575> };



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