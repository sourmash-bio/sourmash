//
// Python 2/3 compatibility: PyInt and PyLong
//

#if (PY_MAJOR_VERSION >= 3)
#define PyInt_Check(arg) PyLong_Check(arg)
#define PyInt_AsLong(arg) PyLong_AsLong(arg)
#define PyInt_FromLong(arg) PyLong_FromLong(arg)
#endif

//
// Python 2/3 compatibility: PyBytes and PyString
// https://docs.python.org/2/howto/cporting.html#str-unicode-unification
//

#include "bytesobject.h"

//
// Python 2/3 compatibility: Module initialization
// http://python3porting.com/cextensions.html#module-initialization
//

#if PY_MAJOR_VERSION >= 3
#define MOD_ERROR_VAL NULL
#define MOD_SUCCESS_VAL(val) val
#define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
#define MOD_DEF(ob, name, doc, methods) \
          static struct PyModuleDef moduledef = { \
            PyModuleDef_HEAD_INIT, name, doc, -1, methods, }; \
          ob = PyModule_Create(&moduledef);
#else
#define MOD_ERROR_VAL
#define MOD_SUCCESS_VAL(val)
#define MOD_INIT(name) void init##name(void)
#define MOD_DEF(ob, name, doc, methods) \
          ob = Py_InitModule3(name, methods, doc);
#endif

//
// Function necessary for Python loading:
//

extern "C" {
    MOD_INIT(_sketch);
}

// Must be first.
#include <Python.h>

#include <string>
#include <set>
#include <exception>
#include <iostream>

#include "third-party/smhasher/MurmurHash3.h"

typedef unsigned long long HashIntoType;
typedef std::set<HashIntoType> CMinHashType;
int _hash_murmur(const std::string& kmer);

typedef struct {
  PyObject_HEAD
  CMinHashType * mins;
  unsigned int num;
  unsigned int ksize;
  long int prime;
  bool is_protein;
} sketch_MinHash_Object;

static
void
sketch_MinHash_dealloc(sketch_MinHash_Object * obj)
{
  delete obj->mins;
  obj->mins = NULL;
  Py_TYPE(obj)->tp_free((PyObject*)obj);
}

static
PyObject *
minhash_add_sequence(sketch_MinHash_Object * me, PyObject * args)
{
  const char * sequence = NULL;
  if (!PyArg_ParseTuple(args, "s", &sequence)) {
    return NULL;
  }
  CMinHashType * mins = me->mins;
  CMinHashType::iterator mins_end;
  
  long int h = 0;

  if (!me->is_protein) {
    std::string seq = sequence;
    for (unsigned int i = 0; i < seq.length() - me->ksize + 1; i++) {
      std::string kmer = seq.substr(i, me->ksize);

      h = _hash_murmur(kmer);
      h = ((h % me->prime) + me->prime) % me->prime;

      // std::cout << "xx h is " << _hash_murmur(kmer) << " for " << kmer << "\n";
      
      // std::cout << "inserting: " << h << " " << me->prime << "\n";

      if (mins->size() == me->num) {
        mins_end = mins->end();
        mins_end--;
        if (h < *mins_end) {
          mins->erase(mins_end);
          mins->insert(h);
        }
      }
      else {
        mins->insert(h);
      }
    }
  }
    
  Py_INCREF(Py_None);
  return Py_None;
}

static
PyObject *
minhash_add_hash(sketch_MinHash_Object * me, PyObject * args)
{
  long int hh;
  if (!PyArg_ParseTuple(args, "l", &hh)) {
    return NULL;
  }
  CMinHashType * mins = me->mins;

  hh = ((hh % me->prime) + me->prime) % me->prime;

  // std::cout << "inserting: " << hh << " " << me->prime << "\n";

  mins->insert(hh);

  if (mins->size() > me->num) {
    CMinHashType::iterator mi = mins->end();
    mi--;
    mins->erase(mi);
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static
PyObject *
minhash_get_mins(sketch_MinHash_Object * me, PyObject * args)
{
  if (!PyArg_ParseTuple(args, "")) {
    return NULL;
  }

  CMinHashType * mins = me->mins;
  PyObject * mins_o = PyList_New(mins->size());

  unsigned int j = 0;
  for (CMinHashType::iterator i = mins->begin(); i != mins->end(); ++i) {
    PyList_SET_ITEM(mins_o, j, PyLong_FromUnsignedLongLong(*i));
    j++;
  }
  return(mins_o);
}

static PyMethodDef sketch_MinHash_methods [] = {
  { "add_sequence",
    (PyCFunction)minhash_add_sequence, METH_VARARGS,
    "Add kmer into MinHash"
  },
  { "add_hash",
    (PyCFunction)minhash_add_hash, METH_VARARGS,
    "Add kmer into MinHash"
  },
  { "get_mins",
    (PyCFunction)minhash_get_mins, METH_VARARGS,
    "Get MinHash signature"
  },
  { NULL, NULL, 0, NULL } // sentinel
};

static
PyObject *
sketch_MinHash_new(PyTypeObject * subtype, PyObject * args, PyObject * kwds)
{
    PyObject * self     = subtype->tp_alloc( subtype, 1 );
    if (self == NULL) {
        return NULL;
    }

    unsigned int _n, _ksize;
    long int _p;
    PyObject * is_protein_o;
    if (!PyArg_ParseTuple(args, "IIlO", &_n, &_ksize, &_p, &is_protein_o)){
      return NULL;
    }
    
    sketch_MinHash_Object * myself = (sketch_MinHash_Object *)self;

    myself->mins = new CMinHashType;
    myself->num = _n;
    myself->ksize = _ksize;
    myself->prime = _p;
    myself->is_protein = PyObject_IsTrue(is_protein_o);

    return self;
}

static PyTypeObject sketch_MinHash_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)        /* init & ob_size */
    "_sketch.MinHash",                    /* tp_name */
    sizeof(sketch_MinHash_Object),         /* tp_basicsize */
    0,                                    /* tp_itemsize */
    (destructor)sketch_MinHash_dealloc,    /* tp_dealloc */
    0,                                    /* tp_print */
    0,                                    /* tp_getattr */
    0,                                    /* tp_setattr */
    0,                                    /* tp_compare */
    0,                                    /* tp_repr */
    0,                                    /* tp_as_number */
    0,                                    /* tp_as_sequence */
    0,                                    /* tp_as_mapping */
    0,                                    /* tp_hash */
    0,                                    /* tp_call */
    0,                                    /* tp_str */
    0,                                    /* tp_getattro */
    0,                                    /* tp_setattro */
    0,                                    /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                   /* tp_flags */
    "A MinHash sketch.",                  /* tp_doc */
    0,                                    /* tp_traverse */
    0,                                    /* tp_clear */
    0,                                    /* tp_richcompare */
    0,                                    /* tp_weaklistoffset */
    0,                                    /* tp_iter */
    0,                                    /* tp_iternext */
    sketch_MinHash_methods,               /* tp_methods */
    0,                                    /* tp_members */
    0,                                    /* tp_getset */
    0,                                         /* tp_base */
    0,                                         /* tp_dict */
    0,                                         /* tp_descr_get */
    0,                                         /* tp_descr_set */
    0,                                         /* tp_dictoffset */
    0,                                         /* tp_init */
    0,                                         /* tp_alloc */
    sketch_MinHash_new,                        /* tp_new */
};

std::string _revcomp(const std::string& kmer)
{
    std::string out = kmer;
    size_t ksize = out.size();

    for (size_t i=0; i < ksize; ++i) {
        char complement;

        switch(kmer[i]) {
        case 'A':
            complement = 'T';
            break;
        case 'C':
            complement = 'G';
            break;
        case 'G':
            complement = 'C';
            break;
        case 'T':
            complement = 'A';
            break;
        default:
            throw std::exception();
            break;
        }
        out[ksize - i - 1] = complement;
    }
    return out;
}

int _hash_murmur(const std::string& kmer)
{
  int out[2];
  HashIntoType h, r;
  uint32_t seed = 0;
  MurmurHash3_x86_32((void *)kmer.c_str(), kmer.size(), seed, &out);
  return out[0];
}

static PyObject * murmur3_forward_hash(PyObject * self, PyObject * args)
{
    const char * kmer;

    if (!PyArg_ParseTuple(args, "s", &kmer)) {
        return NULL;
    }

    return PyLong_FromUnsignedLongLong(_hash_murmur(kmer));
}

static PyMethodDef SketchMethods[] = {
    {
        "hash_murmur3",
        murmur3_forward_hash,
        METH_VARARGS,
        "Calculate the hash value of a k-mer using MurmurHash3 "
        "(with reverse complement)",
    },
    { NULL, NULL, 0, NULL } // sentinel
};

MOD_INIT(_sketch)
{
    if (PyType_Ready( &sketch_MinHash_Type ) < 0) {
        return MOD_ERROR_VAL;
    }

    PyObject * m;

    MOD_DEF(m, "_sketch",
            "interface for the sourmash module low-level extensions",
            SketchMethods);

    if (m == NULL) {
        return MOD_ERROR_VAL;
    }

    Py_INCREF(&sketch_MinHash_Type);
    if (PyModule_AddObject( m, "MinHash",
                            (PyObject *)&sketch_MinHash_Type ) < 0) {
        return MOD_ERROR_VAL;
    }
    return MOD_SUCCESS_VAL(m);
}
