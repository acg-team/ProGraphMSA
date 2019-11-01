#include <cmath>

#ifdef __SSE2__
#include <emmintrin.h>
#endif

template<typename T>
static inline T ls_log(T x)
{
   const T a = 2.44247459618085927548717403238913328776812604856113966238617812902399112761292613763080658235564;
   const T b = -4.2040783745848554315883301529007786406310628696382695994938550046831869207082846248658671;
   const T c = -0.72123729809042963774358701619456664388406302428056983119308906451199556380646306;

   int e;
   T d;

   d = std::frexp(x,&e);
   return a + b/(d-c) + e;
}

#ifdef __SSE2__
static inline void ls_log_add(float *data,float add,unsigned int len)
{
	unsigned int alen = len & ~3u;

   __m128i vx1 = _mm_set1_epi32(126);
   __m128i vx2 = _mm_set1_epi32(0x3f000000);
   __m128 vadd = _mm_set1_ps(add);

   const float a = 2.44247459618085927548717403238913328776812604856113966238617812902399112761292613763080658235564;
   const float b = -4.2040783745848554315883301529007786406310628696382695994938550046831869207082846248658671;
   const float c = -0.72123729809042963774358701619456664388406302428056983119308906451199556380646306;

   __m128 va = _mm_set1_ps(a);
   __m128 vb = _mm_set1_ps(b);
   __m128 vc = _mm_set1_ps(c);

   for(unsigned int i=0; i<alen; i+=4) {
      __m128i vdi = _mm_load_si128((__m128i*)(data+i));
      __m128i vei = _mm_srli_epi32(vdi,23);
      vei = _mm_sub_epi32(vei,vx1);
      __m128 ve = _mm_cvtepi32_ps(vei);

      vdi = _mm_srli_epi32(_mm_slli_epi32(vdi,9),9);
      vdi = _mm_xor_si128(vdi,vx2);

      __m128 vd = reinterpret_cast<__m128>(vdi);

      vd = _mm_div_ps(vb,_mm_sub_ps(vd,vc));
      ve = _mm_add_ps(ve,va);
      vd = _mm_add_ps(vd,vadd);
      vd = _mm_add_ps(vd,ve);

      _mm_store_ps(data+i,vd);
   }
   for(unsigned int i=alen; i<len; ++i) {
      data[i] = ls_log(data[i]) + add;
   }
}

static inline void ls_log_add(double *data,double add,unsigned int len)
{
	unsigned int alen = len & ~1u;

   __m128i vx1 = _mm_set1_epi64((__m64)1022ll);
   __m128i vx2 = _mm_set1_epi64((__m64)0x3fell);
   __m128d vadd = _mm_set1_pd(add);

   const double a = 2.44247459618085927548717403238913328776812604856113966238617812902399112761292613763080658235564;
   const double b = -4.2040783745848554315883301529007786406310628696382695994938550046831869207082846248658671;
   const double c = -0.72123729809042963774358701619456664388406302428056983119308906451199556380646306;

   __m128d va = _mm_set1_pd(a);
   __m128d vb = _mm_set1_pd(b);
   __m128d vc = _mm_set1_pd(c);

   for(unsigned int i=0; i<alen; i+=2) {
      __m128i vdi = _mm_load_si128((__m128i*)(data+i));
      __m128i vei = _mm_srli_epi64(vdi,52);
      vei = _mm_sub_epi32(vei,vx1);
      vei = _mm_xor_si128(vei,_mm_srli_si128(vei,4));
      __m128d ve = _mm_cvtepi32_pd(vei);

      vdi = _mm_srli_epi64(_mm_slli_epi64(vdi,12),12);
      vdi = _mm_xor_si128(vdi,_mm_slli_epi64(vx2,52));

      __m128d vd = reinterpret_cast<__m128d>(vdi);

      vd = _mm_div_pd(vb,_mm_sub_pd(vd,vc));
      ve = _mm_add_pd(ve,va);
      vd = _mm_add_pd(vd,vadd);
      vd = _mm_add_pd(vd,ve);

      _mm_store_pd(data+i,vd);
   }
   for(unsigned int i=alen; i<len; ++i) {
      data[i] = ls_log(data[i]) + add;
   }
}
#else
template<typename T>
static inline void ls_log_add(T *data,T add,unsigned int len)
{
   for(unsigned int i=len; i<len; ++i) {
	  data[i] = ls_log(data[i]) + add;
   }
}
#endif

