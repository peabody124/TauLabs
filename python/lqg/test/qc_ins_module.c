#include <Python.h>
#include "math.h"

#define NPY_NO_DEPRECATED_API 7
#include "numpy/arrayobject.h"
#include "numpy/ndarraytypes.h"

#include <pios_heap.h>
#include <qc_ins.h>

//! State variable
static uintptr_t qcins_handle;

int not_doublevector(PyArrayObject *vec)
{
	if (PyArray_TYPE(vec) != NPY_DOUBLE) {
		PyErr_SetString(PyExc_ValueError,
              "Vector is not a float or double vector.");
		return 1;
	}
	if (PyArray_NDIM(vec) != 1)  {
		PyErr_Format(PyExc_ValueError,
              "Vector is not a 1 dimensional vector (%d).", PyArray_NDIM(vec));
		return 1;  
	}
	return 0;
}

/**
 * parseFloatVecN(vec_in, vec_out, N)
 *
 * @param[in] vec_in the python array to extract elements from
 * @param[out] vec_out float array of the numbers
 * @return true if successful, false if not
 *
 * Convert a python array type to a 3 element float
 * vector.
 */
static bool parseFloatVecN(PyArrayObject *vec_in, float *vec_out, int N)
{
	/* verify it is a valid vector */
	if (not_doublevector(vec_in))
		return false;

	if (PyArray_DIM(vec_in,0) != N) {
		PyErr_Format(PyExc_ValueError, "Vector length incorrect. Received %ld and expected %d", PyArray_DIM(vec_in,0), N);
		return false;
	}

	NpyIter *iter;
	NpyIter_IterNextFunc *iternext;

	/*  create the iterators */
	iter = NpyIter_New(vec_in, NPY_ITER_READONLY, NPY_KEEPORDER,
							 NPY_NO_CASTING, NULL);
	if (iter == NULL)
		goto fail;

	iternext = NpyIter_GetIterNext(iter, NULL);
	if (iternext == NULL) {
		NpyIter_Deallocate(iter);
		goto fail;
	}

	double ** dataptr = (double **) NpyIter_GetDataPtrArray(iter);

	/*  iterate over the arrays */
	int i = 0;
	do {
		vec_out[i++] = **dataptr;
	} while(iternext(iter) && (i < N));

	NpyIter_Deallocate(iter);

	return true;

fail:
	fprintf(stderr, "Parse fail\r\n");
	return false;
}

/**
 * parseFloatVec3(vec_in, vec_out)
 *
 * @param[in] vec_in the python array to extract elements from
 * @param[out] vec_out float array of the numbers
 * @return true if successful, false if not
 *
 * Convert a python array type to a 3 element float
 * vector.
 */
static bool parseFloatVec3(PyArrayObject *vec_in, float *vec_out)
{
  return parseFloatVecN(vec_in, vec_out, 3);
}

/**
 * pack_state put the state information into an array
 */
static PyObject*
pack_state(PyObject* self)
{
	const int N = 15;
	const int nd = 1;
	npy_intp dims[nd] = {N};
	PyArrayObject *state;
	state = (PyArrayObject*) PyArray_SimpleNew(nd, dims, NPY_DOUBLE);
	double *s = (double *) PyArray_DATA(state);

	float p[1];
	float v[3];
	float q[4];
	float rate[3];
	float torque[4];

	bool success = true;
	success &= qcins_get_altitude(qcins_handle, p);
	success &= qcins_get_velocity(qcins_handle, v);
	success &= qcins_get_attitude(qcins_handle, q);
	success &= qcins_get_rate(qcins_handle, rate);
	success &= qcins_get_torque(qcins_handle, torque);

	if (!success)
		return Py_None;

	s[0] = p[0];
	s[1] = v[0];
	s[2] = v[1];
	s[3] = v[2];
	s[4] = q[0];
	s[5] = q[1];
	s[6] = q[2];
	s[7] = q[3];
	s[8] = rate[0];
	s[9] = rate[1];
	s[10] = rate[2];
	s[11] = torque[0];
	s[12] = torque[1];
	s[13] = torque[2];
	s[14] = torque[3];

	return Py_BuildValue("O", state);
}

static PyObject*
init(PyObject* self, PyObject* args)
{
	if (qcins_handle == 0)
		qcins_alloc(&qcins_handle);
	qcins_init(qcins_handle);
	return pack_state(self);
}

static PyObject*
advance(PyObject* self, PyObject* args)
{
	PyArrayObject *vec_control;
	float control[4];
	float dT;

	if (!PyArg_ParseTuple(args, "O!f", &PyArray_Type, &vec_control, &dT))  return NULL;
	if (NULL == vec_control)  return NULL;

	if (!parseFloatVecN(vec_control, control, 4))
		return NULL;

	qcins_predict(qcins_handle, control[0], control[1], control[2], control[3], dT);

	return pack_state(self);
}

static PyObject*
correct(PyObject* self, PyObject* args)
{
	PyArrayObject *vec_gyros;
	PyArrayObject *vec_accels;
	float gyros[3];
	float accels[3];

	if (!PyArg_ParseTuple(args, "O!O!", &PyArray_Type, &vec_gyros, &PyArray_Type, &vec_accels))  return NULL;
	if (NULL == vec_gyros)  return NULL;
	if (NULL == vec_accels)  return NULL;

	if (!parseFloatVec3(vec_gyros, gyros))
		return NULL;
	if (!parseFloatVec3(vec_accels, accels))
		return NULL;

	qcins_correct_accel_gyro(qcins_handle, accels, gyros);
	
	return pack_state(self);
}

static PyMethodDef QcInsMethods[] =
{
	{"init", init, METH_VARARGS, "Reset QC INS state."},
	{"advance", advance, METH_VARARGS, "advance(u, dT) - Advance state 1 time step."},
	{"correct", correct, METH_VARARGS, "correct(g, a) - Correct based on gyro and accel data."},
	{NULL, NULL, 0, NULL}
};
 
PyMODINIT_FUNC
initqc_ins(void)
{
	(void) Py_InitModule("qc_ins", QcInsMethods);
	import_array();
	init(NULL, NULL);
}
