
#pragma once

// #include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
// #include <string.h>

// namespace spiceweasel {

template<typename T, int WINDOW = 3>
class MedianFilter {
public:
  static_assert(WINDOW >= 3, "MedianFilter::window size must be >= 3");
  static_assert(WINDOW % 2, "MedianFilter::window size must be odd");

  MedianFilter() {
    memset(buffer,0,WINDOW*sizeof(T));
    memset(sorted,0,WINDOW*sizeof(T));
  }

  T filter(const T &sample) {
    p = (p + 1) % WINDOW;
    buffer[p] = sample;

    memcpy(sorted, buffer, sizeof(buffer));
    qsort(&sorted, WINDOW, sizeof(T), cmp);

    return sorted[WINDOW / 2];
  }

protected:
  static int cmp(const void *a, const void *b) {
    // return (*(T *)a >= *(T *)b) ? 1 : -1;
    return (*static_cast<T*>(a) >= *static_cast<T*>(b)) ? 1 : -1;
  }

  T buffer[WINDOW]{0.0};
  T sorted[WINDOW]{0.0};
  uint8_t p{0};
};

// } // namespace spiceweasel
