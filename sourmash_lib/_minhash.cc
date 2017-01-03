#include <string>
#include <set>
#include <map>
#include <exception>
#include <iostream>

typedef unsigned long long HashIntoType;
typedef std::set<HashIntoType> CMinHashType;
uint64_t _hash_murmur(const std::string& kmer,
		      const uint32_t seed);

#include "_minhash.hh"
#include "../third-party/smhasher/MurmurHash3.h"

//
// Python 2/3 compatibility: PyInt and PyLong
//

#if (PY_MAJOR_VERSION >= 3)
#define PyInt_Check(arg) PyLong_Check(arg)
#define PyInt_AsLong(arg) PyLong_AsLong(arg)
#define PyInt_FromLong(arg) PyLong_FromLong(arg)
#define Py_TPFLAGS_HAVE_ITER 0
#endif

extern "C" {
    MOD_INIT(_minhash);
}


////

static int _MinHash_len(PyObject *);
static PyObject * _MinHash_concat_inplace(PyObject *, PyObject *);

static PySequenceMethods _MinHash_seqmethods[] = {
    (lenfunc)_MinHash_len, /* sq_length */
    0,      /* sq_concat */
    0,                          /* sq_repeat */
    0,                          /* sq_item */
    0,                          /* sq_slice */
    0,                          /* sq_ass_item */
    0,                          /* sq_ass_slice */
    0,                          /* sq_contains */
    (binaryfunc)_MinHash_concat_inplace, /* sq_inplace_concat */
    0                           /* sq_inplace_repeat */
};


/* Forward declarations */
static void MinHash_dealloc(MinHash_Object * obj);
static PyObject * minhash_add_sequence(MinHash_Object * me, PyObject * args);
static PyObject * minhash_add_protein(MinHash_Object * me, PyObject * args);
static PyObject * minhash_add_hash(MinHash_Object * me, PyObject * args);
static PyObject * minhash_get_mins(MinHash_Object * me, PyObject * args, PyObject * kwargs);
static PyObject * minhash_set_counters(MinHash_Object * me, PyObject * args);
static PyObject * minhash___copy__(MinHash_Object * me, PyObject * args);
static PyObject * minhash_count_common(MinHash_Object * me, PyObject * args);
static PyObject * minhash_compare(MinHash_Object * me, PyObject * args);
static PyObject * minhash_merge(MinHash_Object * me, PyObject * args);
static PyObject * minhash_is_protein(MinHash_Object * me, PyObject * args);
  
static PyMethodDef MinHash_methods [] = {
    {
        "add_sequence",
        (PyCFunction)minhash_add_sequence, METH_VARARGS,
        "Add kmer into MinHash"
    },
    {
        "add_protein",
        (PyCFunction)minhash_add_protein, METH_VARARGS,
        "Add AA kmer into protein MinHash"
    },
    {
        "add_hash",
        (PyCFunction)minhash_add_hash, METH_VARARGS,
        "Add kmer into MinHash"
    },
    {
        "get_mins",
        (PyCFunction)minhash_get_mins, METH_VARARGS | METH_KEYWORDS,
        "Get MinHash signature"
    },
    {
        "set_abundances",
        (PyCFunction)minhash_set_counters, METH_VARARGS,
        "Set abundances for MinHash hashes.",
    },
    {
        "__copy__",
        (PyCFunction)minhash___copy__, METH_VARARGS,
        "Copy this MinHash object",
    },
    {
        "count_common",
        (PyCFunction)minhash_count_common, METH_VARARGS,
        "Get number of hashes in common with other."
    },
    {
        "compare",
        (PyCFunction)minhash_compare, METH_VARARGS,
        "Get the Jaccard similarity between this and other."
    },
    {
        "merge",
        (PyCFunction)minhash_merge, METH_VARARGS,
        "Merge the other MinHash into this one."
    },
    {
        "is_protein",
        (PyCFunction)minhash_is_protein, METH_VARARGS,
        "Return False if a DNA MinHash, True if protein."
    },
    { NULL, NULL, 0, NULL } // sentinel
};


PyDoc_STRVAR(minhash_get_seed_doc,
             "Seed used for hashing (is set when the instance of "
	     "MinHash is created).");
static PyObject * minhash_get_seed(MinHash_Object *self);

/* Getters and setters for MinHash */
static PyGetSetDef MinHash_getseters[] = {
  {"seed", 
   (getter)minhash_get_seed,
   NULL,
   minhash_get_seed_doc,
   NULL},
  {NULL}          /* sentinel */
};


static PyObject * MinHash_new(PyTypeObject * subtype, PyObject * args, PyObject * kwds);


static PyTypeObject MinHash_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)        /* init & ob_size */
    "sourmash_lib._minhash.MinHash",      /* tp_name */
    sizeof(MinHash_Object),               /* tp_basicsize */
    0,                                    /* tp_itemsize */
    (destructor)MinHash_dealloc,          /* tp_dealloc */
    0,                                    /* tp_print */
    0,                                    /* tp_getattr */
    0,                                    /* tp_setattr */
    0,                                    /* tp_compare */
    0,                                    /* tp_repr */
    0,                                    /* tp_as_number */
    _MinHash_seqmethods,                  /* tp_as_sequence */
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
    MinHash_methods,                      /* tp_methods */
    0,                                    /* tp_members */
    MinHash_getseters,                    /* tp_getset */
    0,                                    /* tp_base */
    0,                                    /* tp_dict */
    0,                                    /* tp_descr_get */
    0,                                    /* tp_descr_set */
    0,                                    /* tp_dictoffset */
    0,                                    /* tp_init */
    0,                                    /* tp_alloc */
    MinHash_new,                          /* tp_new */
    0,                                    /* tp_free */
    0                                     /* tp_is_gc */
};




bool check_IsMinHash(PyObject * mh);

PyObject * build_MinHash_Object(KmerMinHash * mh, bool track_abundance = false)
{
    MinHash_Object * obj = (MinHash_Object *) \
                           PyObject_New(MinHash_Object, &MinHash_Type);
    obj->mh = mh;
    obj->track_abundance = track_abundance;

    return (PyObject *) obj;
}

////

static
void
MinHash_dealloc(MinHash_Object * obj)
{
    delete obj->mh;
    obj->mh = NULL;
    Py_TYPE(obj)->tp_free((PyObject*)obj);
}


static
PyObject *
minhash_get_seed(MinHash_Object *self)
{
  return PyLong_FromLong((long)(self->mh->seed));
}

static
PyObject *
minhash_add_sequence(MinHash_Object * me, PyObject * args)
{
    const char * sequence = NULL;
    PyObject * force_o = NULL;
    if (!PyArg_ParseTuple(args, "s|O", &sequence, &force_o)) {
        return NULL;
    }
    KmerMinHash * mh = me->mh;
    bool force = false;
    if (force_o && PyObject_IsTrue(force_o)) {
        force = true;
    }

    try {
        mh->add_sequence(sequence, force);
    } catch (minhash_exception &e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static
PyObject *
minhash_add_protein(MinHash_Object * me, PyObject * args)
{
    const char * sequence = NULL;
    if (!PyArg_ParseTuple(args, "s", &sequence)) {
        return NULL;
    }
    KmerMinHash * mh = me->mh;

    unsigned int ksize = mh->ksize / 3;

    if(strlen(sequence) < ksize) {
        Py_INCREF(Py_None);
        return Py_None;
    }


    if (!mh->is_protein) {
        PyErr_SetString(PyExc_ValueError,
                        "cannot add amino acid sequence to DNA MinHash!");
        return NULL;
    } else {                      // protein
        std::string seq = sequence;
        for (unsigned int i = 0; i < seq.length() - ksize + 1; i ++) {
            std::string aa = seq.substr(i, ksize);

            mh->add_word(aa);
        }
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static
PyObject *
minhash_add_hash(MinHash_Object * me, PyObject * args)
{
    HashIntoType hh;
    if (!PyArg_ParseTuple(args, "K", &hh)) {
        return NULL;
    }

    me->mh->add_hash(hh);

    Py_INCREF(Py_None);
    return Py_None;
}

static
PyObject *
minhash_set_counters(MinHash_Object * me, PyObject * args)
{
    PyObject* dict;
    if (!PyArg_ParseTuple(args, "O!", &PyDict_Type, &dict)) {
        return NULL;
    }

    KmerMinAbundance *mh = (KmerMinAbundance*)me->mh;
    PyObject *key, *value;
    Py_ssize_t pos = 0;

    while (PyDict_Next(dict, &pos, &key, &value)) {
        HashIntoType keyc = 0;
        if (PyLong_Check(key)) {
            keyc = PyLong_AsUnsignedLongLong(key);
        } else if (PyInt_Check(key)) {
            keyc = PyInt_AsLong(key);
        }

        uint64_t valuec = 0;
        if (PyLong_Check(value)) {
            valuec = PyLong_AsUnsignedLongLong(value);
        } else if (PyInt_Check(value)) {
            valuec = PyInt_AsLong(value);
        }

        mh->mins[keyc] = valuec;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static
PyObject *
minhash_get_mins(MinHash_Object * me, PyObject * args, PyObject * kwargs)
{
    PyObject * with_abundance_o = Py_False;
    static const char* const_kwlist[] = {"with_abundance", NULL};
    static char** kwlist = const_cast<char**>(const_kwlist);

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|O", kwlist,
                                     &with_abundance_o)) {
        return NULL;
    }

    bool with_abundance = false;
    if (PyObject_IsTrue(with_abundance_o)) {
        with_abundance = true;
    }

    PyObject * mins_o = NULL;
    if (with_abundance and me->track_abundance) {
        KmerMinAbundance * mh = (KmerMinAbundance*)(me->mh);
        mins_o = PyDict_New();
        for (auto i: mh->mins) {
            PyDict_SetItem(mins_o,
               PyLong_FromUnsignedLongLong(i.first),
               PyLong_FromUnsignedLongLong(i.second));
        }
    } else if (me->track_abundance) {
        KmerMinAbundance * mh = (KmerMinAbundance*)(me->mh);
        mins_o = PyList_New(mh->mins.size());
        unsigned int j = 0;
        for (auto i: mh->mins) {
            PyList_SET_ITEM(mins_o, j, PyLong_FromUnsignedLongLong(i.first));
            j++;
        }
    } else {
        unsigned int j = 0;
        KmerMinHash * mh = me->mh;
        mins_o = PyList_New(mh->mins.size());
        for (auto i: mh->mins) {
            PyList_SET_ITEM(mins_o, j, PyLong_FromUnsignedLongLong(i));
            j++;
        }
    }
    return(mins_o);
}

static int _MinHash_len(PyObject * me)
{
    KmerMinHash * mh = ((MinHash_Object *)me)->mh;
    return mh->num;
}

static PyObject * _MinHash_concat_inplace(PyObject * me_obj,
                                          PyObject * other_obj)
{
    MinHash_Object * me, * other;
    me = (MinHash_Object *) me_obj;
    other = (MinHash_Object *) other_obj;

    try {
        if (me->track_abundance == true) {
            ((KmerMinAbundance*)me->mh)->merge(*(KmerMinAbundance*)other->mh);
        } else {
            me->mh->merge(*other->mh);
        }
    } catch (minhash_exception &e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }

    Py_INCREF(me);
    return (PyObject *) me;
}

static PyObject * minhash___copy__(MinHash_Object * me, PyObject * args)
{
    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    KmerMinHash * mh = me->mh;
    KmerMinHash * new_mh = NULL;

    if (me->track_abundance == true) {
      new_mh = new KmerMinAbundance(mh->num, mh->ksize,
				    mh->is_protein, mh->seed);
         ((KmerMinAbundance*)new_mh)->merge(*(KmerMinAbundance*)mh);
    } else {
      new_mh = new KmerMinHash(mh->num, mh->ksize,
			       mh->is_protein, mh->seed);
         new_mh->merge(*mh);
    }

    return build_MinHash_Object(new_mh, me->track_abundance);
}

static PyObject * minhash_merge(MinHash_Object * me, PyObject * args)
{
    PyObject * other_mh;
    if (!PyArg_ParseTuple(args, "O", &other_mh)) {
        return NULL;
    }
    if (!check_IsMinHash(other_mh)) {
        return NULL;
    }

    KmerMinHash * mh = me->mh;
    KmerMinHash * other = ((MinHash_Object *) other_mh)->mh;

    try {
        if (me->track_abundance == true) {
            ((KmerMinAbundance*)mh)->merge(*(KmerMinAbundance*)other);
        } else {
            mh->merge(*other);
        }
    } catch (minhash_exception &e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }

    Py_INCREF(me);
    return (PyObject *) me;
}

static PyObject * minhash_count_common(MinHash_Object * me, PyObject * args)
{
    PyObject * other_mh;

    if (!PyArg_ParseTuple(args, "O", &other_mh)) {
        return NULL;
    }

    if (!check_IsMinHash(other_mh)) {
        return NULL;
    }
    MinHash_Object * other = (MinHash_Object*) other_mh;

    unsigned int n;
    try {
        if (me->track_abundance) {
            KmerMinAbundance *mh = dynamic_cast<KmerMinAbundance*>(me->mh);
            if (other->track_abundance) {
                KmerMinAbundance *other_mh = dynamic_cast<KmerMinAbundance*>(other->mh);
                n = mh->count_common(*other_mh);
            } else {
                n = mh->count_common(*other->mh);
            }
        } else {
            if (other->track_abundance) {
                KmerMinAbundance *other_mh = dynamic_cast<KmerMinAbundance*>(other->mh);
                n = other_mh->count_common(*me->mh);
            } else {
                n = me->mh->count_common(*other->mh);
            }
        }
    } catch (minhash_exception &e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }
    return PyInt_FromLong(n);
}

static PyObject * minhash_compare(MinHash_Object * me, PyObject * args)
{
    PyObject * other_mh;

    if (!PyArg_ParseTuple(args, "O", &other_mh)) {
        return NULL;
    }

    if (!check_IsMinHash(other_mh)) {
        return NULL;
    }
    MinHash_Object * other = (MinHash_Object*) other_mh;

    unsigned int n;
    unsigned int size;

    try {
        if (me->track_abundance) {
            KmerMinAbundance *mh = dynamic_cast<KmerMinAbundance*>(me->mh);
            if (other->track_abundance) {
                KmerMinAbundance *other_mh = dynamic_cast<KmerMinAbundance*>(other->mh);
                n = mh->count_common(*other_mh);
            } else {
                n = mh->count_common(*other->mh);
            }
            size = mh->mins.size();
        } else {
            if (other->track_abundance) {
                KmerMinAbundance *other_mh = dynamic_cast<KmerMinAbundance*>(other->mh);
                n = other_mh->count_common(*me->mh);
            } else {
                n = me->mh->count_common(*other->mh);
            }
            size = me->mh->mins.size();
        }
    } catch (minhash_exception &e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }

    return PyFloat_FromDouble(float(n) / float(size));
}

static PyObject * minhash_is_protein(MinHash_Object * me, PyObject * args)
{
    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    PyObject * r = NULL;
    if (me->mh->is_protein) {
        r = Py_True;
    } else {
        r = Py_False;
    }
    Py_INCREF(r);

    return r;
}


static
PyObject *
MinHash_new(PyTypeObject * subtype, PyObject * args, PyObject * kwds)
{
    PyObject * self     = subtype->tp_alloc( subtype, 1 );
    if (self == NULL) {
        return NULL;
    }

    unsigned int _n, _ksize;
    uint32_t seed = MINHASH_DEFAULT_SEED;
    PyObject * track_abundance_o = Py_False;
    PyObject * is_protein_o = NULL;
 
    static const char* const_kwlist[] = {"n", "ksize", "is_protein", "track_abundance", "seed", NULL};
    static char** kwlist = const_cast<char**>(const_kwlist);

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "II|OOI", kwlist,
				     &_n, &_ksize, &is_protein_o,
				     &track_abundance_o,
				     &seed)) {
        return NULL;
    }

    MinHash_Object * myself = (MinHash_Object *)self;
    bool is_protein = false;
    if (is_protein_o && PyObject_IsTrue(is_protein_o)) {
        is_protein = true;
    }

    if (PyObject_IsTrue(track_abundance_o)) {
      myself->mh = new KmerMinAbundance(_n, _ksize, is_protein, seed);
        myself->track_abundance = true;
    } else {
      myself->mh = new KmerMinHash(_n, _ksize, is_protein, seed);
        myself->track_abundance = false;
    }

    return self;
}


bool check_IsMinHash(PyObject * mh)
{
    if (!PyObject_TypeCheck(mh, &MinHash_Type)) {
        return false;
    }
    return true;
}


PyDoc_STRVAR(hash_murmur_doc,
             "hash_murmur(string [, seed])\n\n"
             "Compute a hash for a string, optionally using a seed (an integer). "
             "The current default seed is returned by hash_seed().");

static PyObject * hash_murmur(PyObject * self, PyObject * args)
{
    const char * kmer;
    uint32_t seed = MINHASH_DEFAULT_SEED;

      if (!PyArg_ParseTuple(args, "s|I", &kmer, &seed)) {
        return NULL;
    }

    return PyLong_FromUnsignedLongLong(_hash_murmur(kmer, seed));
}

PyDoc_STRVAR(hash_seed_doc,
             "hash_seed() -> int\n\n"
	     "Return the default seed when hashing.");

static PyObject * hash_seed(void *)
{
  return PyLong_FromUnsignedLong((long)MINHASH_DEFAULT_SEED);
}

static PyMethodDef MinHashModuleMethods[] = {
    {
      "hash_murmur", (PyCFunction)hash_murmur,
        METH_VARARGS, hash_murmur_doc,
    },
    {
      "hash_seed", (PyCFunction)hash_seed,
        METH_NOARGS, hash_seed_doc,
    },
    { NULL, NULL, 0, NULL } // sentinel
};


MOD_INIT(_minhash)
{

    PyObject * m;

    MOD_DEF(m, "_minhash",
            "interface for the sourmash module low-level extensions",
            MinHashModuleMethods);

    if (m == NULL) {
        return MOD_ERROR_VAL;
    }

    if (PyType_Ready( &MinHash_Type ) < 0) {
        return MOD_ERROR_VAL;
    }


    Py_INCREF(&MinHash_Type);
    if (PyModule_AddObject( m, "MinHash",
                            (PyObject *)&MinHash_Type ) < 0) {
        return MOD_ERROR_VAL;
    }
    return MOD_SUCCESS_VAL(m);
}

uint64_t _hash_murmur(const std::string& kmer,
		      const uint32_t seed) {
    uint64_t out[2];
    out[0] = 0; out[1] = 0;
    MurmurHash3_x64_128((void *)kmer.c_str(), kmer.size(), seed, &out);
    return out[0];
}
