#include <string>
#include <set>
#include <map>
#include <exception>
#include <iostream>

typedef unsigned long long HashIntoType;
typedef std::set<HashIntoType> CMinHashType;
uint64_t _hash_murmur(const std::string& kmer);

#include "_minhash.hh"
#include "../third-party/smhasher/MurmurHash3.h"


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

PyTypeObject MinHash_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)        /* init & ob_size */
    "_minhash.MinHash",                   /* tp_name */
    sizeof(MinHash_Object),               /* tp_basicsize */
    0,                                    /* tp_itemsize */
    0,                                    /* tp_dealloc */
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
};

bool check_IsMinHash(PyObject * mh);

PyObject * build_MinHash_Object(KmerMinHash * mh)
{
    MinHash_Object * obj = (MinHash_Object *) \
                           PyObject_New(MinHash_Object, &MinHash_Type);
    obj->mh = mh;

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
minhash_get_mins(MinHash_Object * me, PyObject * args)
{
    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    KmerMinHash * mh = me->mh;
    PyObject * mins_o = PyList_New(mh->mins.size());

    unsigned int j = 0;
    for (CMinHashType::iterator i = mh->mins.begin(); i != mh->mins.end(); ++i) {
        PyList_SET_ITEM(mins_o, j, PyLong_FromUnsignedLongLong(*i));
        j++;
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
        me->mh->merge(*other->mh);
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
    KmerMinHash * new_mh = new KmerMinHash(mh->num, mh->ksize, mh->is_protein);
    new_mh->merge(*mh);

    return build_MinHash_Object(new_mh);
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
        mh->merge(*other);
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
        n = me->mh->count_common(*other->mh);
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
        n = me->mh->count_common(*other->mh);
    } catch (minhash_exception &e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }
    size = me->mh->mins.size();

    return PyFloat_FromDouble(float(n) / float(size));
}

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
        (PyCFunction)minhash_get_mins, METH_VARARGS,
        "Get MinHash signature"
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
    { NULL, NULL, 0, NULL } // sentinel
};

static
PyObject *
MinHash_new(PyTypeObject * subtype, PyObject * args, PyObject * kwds)
{
    PyObject * self     = subtype->tp_alloc( subtype, 1 );
    if (self == NULL) {
        return NULL;
    }

    unsigned int _n, _ksize;
    PyObject * is_protein_o = NULL;
    if (!PyArg_ParseTuple(args, "II|O", &_n, &_ksize, &is_protein_o)) {
        return NULL;
    }

    MinHash_Object * myself = (MinHash_Object *)self;
    bool is_protein = false;
    if (is_protein_o && PyObject_IsTrue(is_protein_o)) {
        is_protein = true;
    }

    myself->mh = new KmerMinHash(_n, _ksize, is_protein);

    return self;
}

bool check_IsMinHash(PyObject * mh)
{
    if (!PyObject_TypeCheck(mh, &MinHash_Type)) {
        return false;
    }
    return true;
}


static PyObject * hash_murmur(PyObject * self, PyObject * args)
{
    const char * kmer;

    if (!PyArg_ParseTuple(args, "s", &kmer)) {
        return NULL;
    }

    return PyLong_FromUnsignedLongLong(_hash_murmur(kmer));
}

static PyMethodDef MinHashModuleMethods[] = {
    {
        "hash_murmur",     hash_murmur,
        METH_VARARGS,       "",
    },
    { NULL, NULL, 0, NULL } // sentinel
};

MOD_INIT(_minhash)
{
    MinHash_Type.tp_methods = MinHash_methods;
    MinHash_Type.tp_dealloc = (destructor)MinHash_dealloc;
    MinHash_Type.tp_new = MinHash_new;

    if (PyType_Ready( &MinHash_Type ) < 0) {
        return MOD_ERROR_VAL;
    }

    PyObject * m;

    MOD_DEF(m, "_minhash",
            "interface for the sourmash module low-level extensions",
            MinHashModuleMethods);

    if (m == NULL) {
        return MOD_ERROR_VAL;
    }

    Py_INCREF(&MinHash_Type);
    if (PyModule_AddObject( m, "MinHash",
                            (PyObject *)&MinHash_Type ) < 0) {
        return MOD_ERROR_VAL;
    }
    return MOD_SUCCESS_VAL(m);
}

uint64_t _hash_murmur(const std::string& kmer) {
    uint64_t out[2];
    out[0] = 0; out[1] = 0;
    uint32_t seed = 42;
    MurmurHash3_x64_128((void *)kmer.c_str(), kmer.size(), seed, &out);
    return out[0];
}
