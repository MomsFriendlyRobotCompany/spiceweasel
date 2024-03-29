

#include "cMathlib.h"

#define BMDEBUG 0

//--- static members ---

static char errMsg[ERROR_STRING_SIZE];
std::vector<cBaseMath *> cBaseMath::cache;
int cBaseMath::keyNum = 10;

/*!
Constructor for cBaseMath.
 */
cBaseMath::cBaseMath(void) {
  key  = 0;
  p    = nullptr;
  name = nullptr;
}

/*!
Base class destructor. It deletes the data array and name, and sets
 their pointers to NULL.
 */
cBaseMath::~cBaseMath(void) {
  if (p) delete[] p;
  if (name) delete[] name;
  p    = nullptr;
  name = nullptr;
}

/*!
This function sets the name of the cMathlib object
 which is used for debugging.
 \param b a string to identify the cMathlib object.
 */
void cBaseMath::setName(const char *b) {
  int size = 0;

  if (name) delete[] name;

  if (b == nullptr) {
    name    = new char[1];
    name[0] = '\0';
  }
  else {
    size = strlen(b);
    name = new char[size + 1];
    if (name == nullptr) {
      snprintf(errMsg, 255, "ERROR: couldn't allocate memory for name %s", b);
      throw cMLError(cMLError::FATAL, errMsg);
    }
    name[size] = '\0';
    strcpy(name, b);
  }
}

/*
 Returns a pointer to a temporary cMathlib object.
 \param r number of rows.
 \param c number of columns.
 \return a temporary cMatrix.
 */
cMatrix *cBaseMath::getTmp(int r, int c) {
  cBaseMath *b       = nullptr;
  cMatrix *m         = nullptr;
  unsigned int index = 0;

  for (index = 0; index < cache.size(); index++) {
    b = cache[index];
    if (b->type == cBaseMath::ML_MATRIX && b->available) {
      m = (cMatrix *)b;
      if ((m->r == r) && (m->c == c)) {
        m->available = false;
        // printf(" using matrix[%d](%d,%d)\n",m->getKey(),r,c);
        return m; // break; // get out of for loop
      }
      else {
        // printf(" wrong size %s %d x %d want %d x
        // %d\n",m->name,m->r,m->c,r,c);
        m = nullptr;
      }
    }
  }

  // no matrix found, create a new one and return it.
  if (m == nullptr) {
    m = new cMatrix(r, c, "temp matrix");
    // m->type == cBaseMath::ML_MATRIX;
    m->setKey(keyNum++);
    m->available = false;
#if BMDEBUG
    printf(" + matrix[%d](%d,%d)\n", m->getKey(), r, c);
#endif
    cache.push_back(m);
  }

  return m;
}

/*!
Returns a pointer to a temporary cMathlib object.
 \param size length of the vector.
 \return a temporary cVector.
 */
cVector *cBaseMath::getTmp(int size) {
  cBaseMath *b       = nullptr;
  cVector *v         = nullptr;
  unsigned int index = 0;

  for (index = 0; index < cache.size(); index++) {
    b = (cBaseMath *)cache[index];
    if (b->type == cBaseMath::ML_VECTOR && b->available) {
      v = (cVector *)b;
      if (v->getSize() == size) {
        v->available = false;
        // printf(" using vector[%d](%d)\n",v->getKey(),size);
        break; // get out of for loop
      }
      else {
        v = nullptr;
      }
    }
  }

  if (v == nullptr) {
    v = new cVector(size, "temp");
    v->setKey(keyNum++);
    v->available = false;
#if BMDEBUG
    printf(" + vector[%d](%d)\n", v->getKey(), size);
#endif
    cache.push_back(v);
  }

  return v;
}

/*!
Returns a pointer to a temporary cMathlib object.
 \return a temporary cQuaternion.
 */
cQuaternion *cBaseMath::getTmp(void) {
  cBaseMath *b       = nullptr;
  cQuaternion *q     = nullptr;
  unsigned int index = 0;

  for (index = 0; index < cache.size(); index++) {
    b = (cBaseMath *)cache[index];
    if (b->type == cBaseMath::ML_QUATERNION && b->available) {
      q            = (cQuaternion *)b;
      q->available = false;
      // printf(" using vector[%d](%d)\n",v->getKey(),size);
      break; // get out of for loop
    }
  }

  if (q == nullptr) {
    q = new cQuaternion("temp");
    q->setKey(keyNum++);
    q->available = false;
#if BMDEBUG
    printf(" + quaternion[%d]\n", q->getKey());
#endif
    cache.push_back(q);
  }

  return q;
}

/*
 Returns a pointer to a temporary cMathlib object.
 \param size length of array.
 \return a temporary cArray.
 */
cArray *cBaseMath::getTmpArray(int size) {
  cBaseMath *b       = nullptr;
  cArray *a          = nullptr;
  unsigned int index = 0;

  for (index = 0; index < cache.size(); index++) {
    b = (cBaseMath *)cache[index];
    if (b->type == cBaseMath::ML_ARRAY && b->available) {
      a = (cArray *)b;
      if (a->getSize() >= size) {
        a->available = false;
        // printf(" using array[%d](%d)\n",a->getKey(),size);
        break; // get out of for loop
      }
      else {
        a = nullptr;
      }
    }
  }

  if (a == nullptr) {
    a = new cArray(size);
    a->setKey(keyNum++);
    a->available = false;
#if BMDEBUG
    printf(" + array[%d](%d)\n", a->getKey(), size);
#endif
    cache.push_back(a);
  }

  memset(a->array, 0, sizeof(int) * size);

  return a;
}

/*!
This function releases a temporary cMathlib object
 so that it is available to be reused later. Temporary
 objects are identified by the presence of a key
 value.
 */
void cBaseMath::releaseTmp(void) const {
  if (key) {
    const_cast<cBaseMath *>(this)->available = true;
  }
}
