/***************************************************************************************************
*
***************************************************************************************************/
#ifndef _EXTRACTOR_KALMAN_H
#define _EXTRACTOR_KALMAN_H

/* define a kalman filter to track a 1d variable */
typedef struct _alg_kf_1d {
  /* the current state var X(k),updated at K moment */
  VSMat     t_X; /* state vector */
  VSMat     t_XP;
  float     f_X, f_XP;

  /* state propagate matrix, current an I matrix,*/
  VSMat     t_A;  /* system transition matrix */
  VSMat     t_AT;
  float     f_A, f_AT;

  VSMat     t_b;
  float     f_b;

  VSMat     t_P; /* cov for X */
  VSMat     t_PP; /* predicted cov for X */
  float     f_P, f_PP;

  VSMat     t_H; /* observing matrix */
  VSMat     t_HT;
  float     f_H, f_HT;


  VSMat     t_R; /* cov for noise in observing equation */
  float     f_R;

  VSMat     t_Q; /* cov for noisy in system transi equa.*/
  VSMat     t_QT;
  float     f_Q, f_QT;

  VSMat     t_I;
  float     f_I;

} ALG_KF_1D;

/* define a kalman filter to track a 2d variable */
typedef struct _alg_kf_2d {
  /* the current state var X(k),updated at K moment */
  VSMat     t_X; /* state vector */
  float     af_X[4];
  VSMat     t_XP;
  float     af_XP[4];

  /* state propagate matrix, current an I matrix,*/
  VSMat     t_A;  /* system transition matrix */
  float     af_A[16];
  VSMat     t_AT;
  float     af_AT[16];

  VSMat     t_b;
  float     af_b[4];

  VSMat     t_P; /* cov for X */
  VSMat     t_PP; /* predicted cov for X */
  float     af_P[16];
  float     af_PP[16];

  VSMat     t_H; /* observing matrix */
  VSMat     t_HT;
  float     af_H[8], af_HT[8];

  VSMat     t_R; /* cov for noise in observing equation */
  float     af_R[4];

  VSMat     t_Q; /* cov for noisy in system transi equa.*/
  VSMat     t_QT;
  float     af_Q[16], af_QT[16];

  VSMat     t_I;
  float     af_I[16];

} ALG_KF_2D;

typedef struct _alg_kf_3d {
  /* the current state var X(k),updated at K moment */
  VSMat     t_X; /* state vector */
  float     af_X[3];
  VSMat     t_XP;
  float     af_XP[3];

  /* state propagate matrix, current an I matrix,*/
  VSMat     t_A;  /* system transition matrix */
  float     af_A[9];
  VSMat     t_AT;
  float     af_AT[9];

  VSMat     t_b;
  float     af_b[3];

  VSMat     t_P; /* cov for X */
  VSMat     t_PP; /* predicted cov for X */
  float     af_P[9];
  float     af_PP[9];

  VSMat     t_H; /* observing matrix */
  VSMat     t_HT;
  float     af_H[9], af_HT[9];

  VSMat     t_R; /* cov for noise in observing equation */
  float     af_R[9];

  VSMat     t_Q; /* cov for noisy in system transi equa.*/
  VSMat     t_QT;
  float     af_Q[9], af_QT[9];

  VSMat     t_I;
  float     af_I[9];

} ALG_KF_3D;

typedef struct _alg_kf_st {

  int         dw_order;

  /* the current state var X(k),updated at K moment */
  VSMat       t_X; /* state vector */
  float      *pf_X;
  VSMat       t_XP;
  float      *pf_XP;

  /* state propagate matrix, current an I matrix,*/
  VSMat       t_A;  /* system transition matrix */
  float      *pf_A;
  VSMat       t_AT;
  float      *pf_AT;

  VSMat       t_P; /* cov for X */
  VSMat       t_PP; /* predicted cov for X */
  float      *pf_P;
  float      *pf_PP;

  VSMat       t_H; /* observing matrix */
  VSMat       t_HT;
  float      *pf_H, *pf_HT;

  VSMat       t_R; /* cov for noise in observing equation */
  float      *pf_R;

  VSMat       t_Q; /* cov for noisy in system transi equa.*/
  VSMat       t_QT;
  float      *pf_Q, *pf_QT;

  VSMat       t_I;
  float      *pf_I;
} ALG_KF_ST;

//! the inner function for ed-line module
/**
*  predict state vectors,covariance matrix of state vectors and measurement vectors
*   Y(n|n-1) = A*Y(n-1|n-1)
*   X(n|n-1) = H*Y(n|n-1)
*   Ryy(n|n-1) = A*Ryy()*A'+Ru;
*
*
* @param[in] pt_out the handle of algorithm module
* @param[in] pt_obs the structure containing edge pixels chains
* @param[out] pt_Xp the orientation(horz&vert) of the starting pixel
*
* @return null
*/
int BuildKalmanFilter_1D(ALG_KF_1D *pt_out);

//! the inner function for ed-line module
/**
*  predict state vectors,covariance matrix of state vectors and measurement vectors
*   Y(n|n-1) = A*Y(n-1|n-1)
*   X(n|n-1) = H*Y(n|n-1)
*   Ryy(n|n-1) = A*Ryy()*A'+Ru;
*
*
* @param[in] pt_out the handle of algorithm module
* @param[in] pt_obs the structure containing edge pixels chains
* @param[out] pt_Xp the orientation(horz&vert) of the starting pixel
*
* @return null
*/
int InitKalmanFilter_1D(ALG_KF_1D *pt_out, float *pf_init_sta, float f_P, float f_PP);

//! the inner function for ed-line module
/**
*  predict state vectors,covariance matrix of state vectors and measurement vectors
*   Y(n|n-1) = A*Y(n-1|n-1)
*   X(n|n-1) = H*Y(n|n-1)
*   Ryy(n|n-1) = A*Ryy()*A'+Ru;
*
*
* @param[in] pt_out the handle of algorithm module
* @param[in] pt_obs the structure containing edge pixels chains
* @param[out] pt_Xp the orientation(horz&vert) of the starting pixel
*
* @return null
*/
int RunKalmanFilter_1D(ALG_KF_1D *pt_out, float *pf_obs, float *pf_pred, float *pf_cur);

//! the inner function for ed-line module
/**
*  predict state vectors,covariance matrix of state vectors and measurement vectors
*   Y(n|n-1) = A*Y(n-1|n-1)
*   X(n|n-1) = H*Y(n|n-1)
*   Ryy(n|n-1) = A*Ryy()*A'+Ru;
*
*
* @param[in] pt_out the handle of algorithm module
* @param[out] pt_Xp the orientation(horz&vert) of the starting pixel
*
* @return null
*/
int RunPredictionKalmanFilter_1D(ALG_KF_1D *pt_out, float *pf_pred);

//! the inner function for ed-line module
/**
*  predict state vectors,covariance matrix of state vectors and measurement vectors
*   Y(n|n-1) = A*Y(n-1|n-1)
*   X(n|n-1) = H*Y(n|n-1)
*   Ryy(n|n-1) = A*Ryy()*A'+Ru;
*
*
* @param[in] pt_out the handle of algorithm module
* @param[in] pt_obs the structure containing edge pixels chains
* @param[out] pt_Xp the orientation(horz&vert) of the starting pixel
*
* @return null
*/
int RunUpdateKalmanFilter_1D(ALG_KF_1D *pt_out, float *pf_obs);

//! the inner function for ed-line module
/**
*  predict state vectors,covariance matrix of state vectors and measurement vectors
*   Y(n|n-1) = A*Y(n-1|n-1)
*   X(n|n-1) = H*Y(n|n-1)
*   Ryy(n|n-1) = A*Ryy()*A'+Ru;
*
*
* @param[in] pt_out the handle of algorithm module
* @param[in] pt_obs the structure containing edge pixels chains
* @param[out] pt_Xp the orientation(horz&vert) of the starting pixel
*
* @return null
*/
int BuildKalmanFilter_2D(ALG_KF_2D *pt_out);

//! the inner function for ed-line module
/**
*  predict state vectors,covariance matrix of state vectors and measurement vectors
*   Y(n|n-1) = A*Y(n-1|n-1)
*   X(n|n-1) = H*Y(n|n-1)
*   Ryy(n|n-1) = A*Ryy()*A'+Ru;
*
*
* @param[in] pt_out the handle of algorithm module
* @param[in] pt_obs the structure containing edge pixels chains
* @param[out] pt_Xp the orientation(horz&vert) of the starting pixel
*
* @return null
*/
int InitKalmanFilter_2D(ALG_KF_2D *pt_out, VSPoint2f *pt_init_state, float *af_P, float *af_PP);

//! the inner function for ed-line module
/**
*  predict state vectors,covariance matrix of state vectors and measurement vectors
*   Y(n|n-1) = A*Y(n-1|n-1)
*   X(n|n-1) = H*Y(n|n-1)
*   Ryy(n|n-1) = A*Ryy()*A'+Ru;
*
*
* @param[in] pt_out the handle of algorithm module
* @param[in] pt_obs the structure containing edge pixels chains
* @param[out] pt_Xp the orientation(horz&vert) of the starting pixel
*
* @return null
*/
int RunKalmanFilter_2D(ALG_KF_2D *pt_out, VSPoint2f *pt_obs, VSPoint2f *pt_pred, VSPoint2f *pt_cur);

//! the inner function for ed-line module
/**
*  predict state vectors,covariance matrix of state vectors and measurement vectors
*   Y(n|n-1) = A*Y(n-1|n-1)
*   X(n|n-1) = H*Y(n|n-1)
*   Ryy(n|n-1) = A*Ryy()*A'+Ru;
*
*
* @param[in] pt_out the handle of algorithm module
* @param[out] pt_Xp the orientation(horz&vert) of the starting pixel
*
* @return null
*/
int RunPredictionKalmanFilter_2D(ALG_KF_2D *pt_out, VSPoint2f *pt_pred);

//! the inner function for ed-line module
/**
*  predict state vectors,covariance matrix of state vectors and measurement vectors
*   Y(n|n-1) = A*Y(n-1|n-1)
*   X(n|n-1) = H*Y(n|n-1)
*   Ryy(n|n-1) = A*Ryy()*A'+Ru;
*
*
* @param[in] pt_out the handle of algorithm module
* @param[in] pt_obs the structure containing edge pixels chains
* @param[out] pt_Xp the orientation(horz&vert) of the starting pixel
*
* @return null
*/
int RunUpdateKalmanFilter_2D(ALG_KF_2D *pt_out, VSPoint2f *pt_Xt);

//! the inner function for ed-line module
/**
*  predict state vectors,covariance matrix of state vectors and measurement vectors
*   Y(n|n-1) = A*Y(n-1|n-1)
*   X(n|n-1) = H*Y(n|n-1)
*   Ryy(n|n-1) = A*Ryy()*A'+Ru;
*
*
* @param[in] pt_out the handle of algorithm module
* @param[in] pt_obs the structure containing edge pixels chains
* @param[out] pt_Xp the orientation(horz&vert) of the starting pixel
*
* @return null
*/
int BuildKalmanFilter_3D(ALG_KF_3D *pt_out);

//! the inner function for ed-line module
/**
*  predict state vectors,covariance matrix of state vectors and measurement vectors
*   Y(n|n-1) = A*Y(n-1|n-1)
*   X(n|n-1) = H*Y(n|n-1)
*   Ryy(n|n-1) = A*Ryy()*A'+Ru;
*
*
* @param[in] pt_out the handle of algorithm module
* @param[in] pt_obs the structure containing edge pixels chains
* @param[out] pt_Xp the orientation(horz&vert) of the starting pixel
*
* @return null
*/
int InitKalmanFilter_3D(ALG_KF_3D *pt_out, VSPoint3f *pt_init_state, float *af_P, float *af_PP);

//! the inner function for ed-line module
/**
*  predict state vectors,covariance matrix of state vectors and measurement vectors
*   Y(n|n-1) = A*Y(n-1|n-1)
*   X(n|n-1) = H*Y(n|n-1)
*   Ryy(n|n-1) = A*Ryy()*A'+Ru;
*
*
* @param[in] pt_out the handle of algorithm module
* @param[in] pt_obs the structure containing edge pixels chains
* @param[out] pt_Xp the orientation(horz&vert) of the starting pixel
*
* @return null
*/
int RunKalmanFilter_3D(ALG_KF_3D *pt_out, VSPoint3f *pt_obs, VSPoint3f *pt_pred, VSPoint3f *pt_cur);

//! the inner function for ed-line module
/**
*  predict state vectors,covariance matrix of state vectors and measurement vectors
*   Y(n|n-1) = A*Y(n-1|n-1)
*   X(n|n-1) = H*Y(n|n-1)
*   Ryy(n|n-1) = A*Ryy()*A'+Ru;
*
*
* @param[in] pt_out the handle of algorithm module
* @param[in] pt_obs the structure containing edge pixels chains
* @param[out] pt_Xp the orientation(horz&vert) of the starting pixel
*
* @return null
*/
int RunPredictionKalmanFilter_3D(ALG_KF_3D *pt_out, VSPoint3f *pt_pred);

//! the inner function for ed-line module
/**
*  predict state vectors,covariance matrix of state vectors and measurement vectors
*   Y(n|n-1) = A*Y(n-1|n-1)
*   X(n|n-1) = H*Y(n|n-1)
*   Ryy(n|n-1) = A*Ryy()*A'+Ru;
*
*
* @param[in] pt_out the handle of algorithm module
* @param[in] pt_obs the structure containing edge pixels chains
* @param[out] pt_Xp the orientation(horz&vert) of the starting pixel
*
* @return null
*/
int RunUpdateKalmanFilter_3D(ALG_KF_3D *pt_out, VSPoint3f *pt_obs);

#endif
