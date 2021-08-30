/***************************************************************************************************
*
*
***************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#if defined(_PC)||defined(_ARM)
#include <memory.h>
#endif
#include <math.h>

#include "base/datatypes.h"
#include "base/defines.h"
#include "math/math.h"

#include "extractor_kalman.h"

#ifdef _DSP_C66X_
#include "bmvs_math_c66x.h"
#endif

int BuildKalmanFilter_1D(ALG_KF_1D *pt_out) {
  VSSetMATHeaderWithType(&pt_out->t_X, 1, 1, 1, sizeof(float), (void *)&pt_out->f_X, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_XP, 1, 1, 1, sizeof(float), (void *)&pt_out->f_XP, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_A, 1, 1, 1, sizeof(float), (void *)&pt_out->f_A, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_AT, 1, 1, 1, sizeof(float), (void *)&pt_out->f_AT, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_b, 1, 1, 1, sizeof(float), (void *)&pt_out->f_b, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_P, 1, 1, 1, sizeof(float), (void *)&pt_out->f_P, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_PP, 1, 1, 1, sizeof(float), (void *)&pt_out->f_PP, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_H, 1, 1, 1, sizeof(float), (void *)&pt_out->f_H, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_HT, 1, 1, 1, sizeof(float), (void *)&pt_out->f_HT, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_R, 1, 1, 1, sizeof(float), (void *)&pt_out->f_R, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_Q, 1, 1, 1, sizeof(float), (void *)&pt_out->f_Q, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_QT, 1, 1, 1, sizeof(float), (void *)&pt_out->f_QT, T_FLOAT);

  pt_out->f_I = 1.f;
  VSSetMATHeaderWithType(&pt_out->t_I, 1, 1, 1, sizeof(float), (void *)&pt_out->f_I, T_FLOAT);
  return VS_SUCCESS;
}

int InitKalmanFilter_1D(ALG_KF_1D *pt_out, float *pf_init_state, float f_P, float f_PP) {
  float f_pred, f_cur;

  /* use current observation as prediction */
  pt_out->f_X = *pf_init_state;
  pt_out->f_XP = *pf_init_state;
  pt_out->f_P = f_P;
  pt_out->f_PP = f_PP;
  RunKalmanFilter_1D(pt_out, pf_init_state, &f_pred, &f_cur);

  return VS_SUCCESS;
}

int RunKalmanFilter_1D(ALG_KF_1D *pt_out, float *pf_obs, float *pf_pred, float *pf_cur) {
  /*
  * for convenience for the whole detect(based on prediction)-match-track(based on current detection) framework,
  * at first,kalman tracking is updated for current frame based on prediction at previous moment and observation at
  * current moment,meanwhile currently updated state as outputing model.
  * then, prediction is performed for usage at last frame
  */
  RunUpdateKalmanFilter_1D(pt_out, pf_obs);

  *pf_cur = pt_out->f_X;

  RunPredictionKalmanFilter_1D(pt_out, pf_pred);

  return VS_SUCCESS;
}

int RunUpdateKalmanFilter_1D(ALG_KF_1D *pt_out, float *pf_obs) {
  VSMat     tPPHT, t_HXP, tKH, tIKH;
  VSMat     tHPP, t_Xt;
  VSMat     tHPPHT, tHPPHTR, tInv;
  VSMat     tK, tInnov, tInnovK;
  float     fPPHT, fHPP, fHPPHT, fHPPHTR, fInv;
  float     fK, fHXP, f_XT, fInnov, fInnovK, fKH, fIKH;

  /* 1 kalman gain: K = PP*H'*inv(H*PP*H'+R); */
  VSSetMATHeaderWithType(&tPPHT, 1, 1, 1, sizeof(float), (void *)&fPPHT, T_FLOAT);
  VSMultMat(&pt_out->t_PP, &pt_out->t_HT, &tPPHT);

  VSSetMATHeaderWithType(&tHPP, 1, 1, 1, sizeof(float), (void *)&fHPP, T_FLOAT);
  VSMultMat(&pt_out->t_H, &pt_out->t_PP, &tHPP);

  VSSetMATHeaderWithType(&tHPPHT, 1, 1, 1, sizeof(float), (void *)&fHPPHT, T_FLOAT);
  VSMultMat(&tHPP, &pt_out->t_HT, &tHPPHT);

  VSSetMATHeaderWithType(&tHPPHTR, 1, 1, 1, sizeof(float), (void *)&fHPPHTR, T_FLOAT);
  VSAddMat(&tHPPHT, &pt_out->t_R, &tHPPHTR);

  VSSetMATHeaderWithType(&tInv, 1, 1, 1, sizeof(float), (void *)&fInv, T_FLOAT);
  ((float *)tInv.data)[0] = 1 / ((float *)tHPPHTR.data)[0];

  VSSetMATHeaderWithType(&tK, 1, 1, 1, sizeof(float), (void *)&fK, T_FLOAT);
  VSMultMat(&tPPHT, &tInv, &tK);

  VSSetMATHeaderWithType(&t_HXP, 1, 1, 1, sizeof(float), (void *)&fHXP, T_FLOAT);
  VSMultMat(&pt_out->t_H, &pt_out->t_XP, &t_HXP);

  f_XT = *pf_obs;
  VSSetMATHeaderWithType(&t_Xt, 1, 1, 1, sizeof(float), (void *)&f_XT, T_FLOAT);
  VSSetMATHeaderWithType(&tInnov, 1, 1, 1, sizeof(float), (void *)&fInnov, T_FLOAT);
  VSSubMat(&t_Xt, &t_HXP, &tInnov);

  /* Kalman gain multipled by innovative  */
  VSSetMATHeaderWithType(&tInnovK, 1, 1, 1, sizeof(float), (void *)&fInnovK, T_FLOAT);
  VSMultMat(&tK, &tInnov, &tInnovK);

  /* update state vector: Xt = xp + K*([cc(t),cr(t)]' - H*xp)*/
  VSAddMat(&tInnovK, &pt_out->t_XP, &pt_out->t_X);

  /*3. P = (eye(4)-K*H)*PP */
  VSSetMATHeaderWithType(&tKH, 1, 1, 1, sizeof(float), (void *)&fKH, T_FLOAT);
  VSMultMat(&tK, &pt_out->t_H, &tKH);

  VSSetMATHeaderWithType(&tIKH, 1, 1, 1, sizeof(float), (void *)&fIKH, T_FLOAT);
  VSSubMat(&pt_out->t_I, &tKH, &tIKH);

  VSMultMat(&tIKH, &pt_out->t_PP, &pt_out->t_P);

  return VS_SUCCESS;
}

int RunPredictionKalmanFilter_1D(ALG_KF_1D *pt_out, float *pf_pred) { //the second is no use
  VSMat         tAX, tAP, tAPAT;
  float         fAX, fAP, fAPAT;
  /*
  *  so the predictions include:
  *    Y(n|n-1) = A*Y(n-1|n-1)
  *    X(n|n-1) = H*Y(n|n-1)
  *    PP(n|n-1) = A*Pyy()*A' + Q';
  */
  /* the predicted states vector: Y(n|n-1) = A*Y(n-1|n-1) */

  VSSetMATHeaderWithType(&tAX, 1, 1, 1, sizeof(float), (void *)&fAX, T_FLOAT);
  VSMultMat(&pt_out->t_A, &pt_out->t_X, &tAX);

  VSAddMat(&tAX, &pt_out->t_b, &pt_out->t_XP);
  *pf_pred = pt_out->f_XP;

  /* Ry(n|n-1) */
  VSSetMATHeaderWithType(&tAP, 1, 1, 1, sizeof(float), (void *)&fAP, T_FLOAT);
  VSMultMat(&pt_out->t_A, &pt_out->t_P, &tAP);

  VSSetMATHeaderWithType(&tAPAT, 1, 1, 1, sizeof(float), (void *)&fAPAT, T_FLOAT);
  VSMultMat(&tAP, &pt_out->t_AT, &tAPAT);

  VSAddMat(&tAPAT, &pt_out->t_QT, &pt_out->t_PP);

  return VS_SUCCESS;
}

int BuildKalmanFilter_2D(ALG_KF_2D *pt_out) {
  VSSetMATHeaderWithType(&pt_out->t_X, 1, 4, 1, sizeof(float), (void *)&pt_out->af_X, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_XP, 1, 4, 1, sizeof(float), (void *)&pt_out->af_XP, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_A, 4, 4, 1, sizeof(float), (void *)&pt_out->af_A, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_AT, 4, 4, 1, sizeof(float), (void *)&pt_out->af_AT, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_b, 1, 1, 1, sizeof(float), (void *)&pt_out->af_b, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_P, 4, 4, 1, sizeof(float), (void *)&pt_out->af_P, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_PP, 4, 4, 1, sizeof(float), (void *)&pt_out->af_PP, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_H, 4, 2, 1, sizeof(float), (void *)&pt_out->af_H, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_HT, 2, 4, 1, sizeof(float), (void *)&pt_out->af_HT, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_R, 2, 2, 1, sizeof(float), (void *)&pt_out->af_R, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_Q, 4, 4, 1, sizeof(float), (void *)&pt_out->af_Q, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_QT, 4, 4, 1, sizeof(float), (void *)&pt_out->af_QT, T_FLOAT);

  pt_out->af_I[0] = 1.f; pt_out->af_I[5] = 1.f; pt_out->af_I[10] = 1.f; pt_out->af_I[15] = 1.f;
  VSSetMATHeaderWithType(&pt_out->t_I, 4, 4, 1, sizeof(float), (void *)&pt_out->af_I, T_FLOAT);

  return VS_SUCCESS;
}

int InitKalmanFilter_2D(ALG_KF_2D *pt_out, VSPoint2f *pt_init_state, float *af_P, float *af_PP) {
  VSPoint2f t_pred, t_cur;

  /* use current observation as prediction */
  pt_out->af_X[0] = pt_init_state->x;
  pt_out->af_X[1] = pt_init_state->y;

  pt_out->af_XP[0] = pt_init_state->x;
  pt_out->af_XP[1] = pt_init_state->y;

  memcpy(pt_out->af_P, af_P, sizeof(pt_out->af_P));
  memcpy(pt_out->af_PP, af_PP, sizeof(pt_out->af_PP));

  RunKalmanFilter_2D(pt_out, pt_init_state, &t_pred, &t_cur);

  return VS_SUCCESS;
}

int RunKalmanFilter_2D(ALG_KF_2D *pt_out, VSPoint2f *pt_obs, VSPoint2f *pt_pred, VSPoint2f *pt_cur) {
  /*
  * for convenience for the whole detect(based on prediction)-match-track(based on current detection) framework,
  * at first,kalman tracking is updated for current frame based on prediction at previous moment and observation at
  * current moment,meanwhile currently updated state as outputing model.
  * then, prediction is performed for usage at last frame
  */
  RunUpdateKalmanFilter_2D(pt_out, pt_obs);

  pt_cur->x = pt_out->af_X[0];
  pt_cur->y = pt_out->af_X[1];

  RunPredictionKalmanFilter_2D(pt_out, pt_pred);

  return VS_SUCCESS;
}

int RunPredictionKalmanFilter_2D(ALG_KF_2D *pt_out, VSPoint2f *pt_pred) { //the second is no use
  VSMat         tAX, tAP, tAPAT;
  float         fAX[4], fAP[4 * 4], fAPAT[4 * 4];
  /*
  *  so the predictions include:
  *    Y(n|n-1) = A*Y(n-1|n-1)
  *    X(n|n-1) = H*Y(n|n-1)
  *    PP(n|n-1) = A*Pyy()*A' + Q';
  */
  /* the predicted states vector: Y(n|n-1) = A*Y(n-1|n-1) */
  VSSetMATHeaderWithType(&tAX, 1, 4, 1, sizeof(float), (void *)fAX, T_FLOAT);
  VSMultMat(&pt_out->t_A, &pt_out->t_X, &tAX);
  VSAddMat(&tAX, &pt_out->t_b, &pt_out->t_XP);
  pt_pred->x = pt_out->af_XP[0];
  pt_pred->y = pt_out->af_XP[1];

  /* Ry(n|n-1) */
  VSSetMATHeaderWithType(&tAP, 4, 4, 1, sizeof(float), (void *)fAP, T_FLOAT);
  VSMultMat(&pt_out->t_A, &pt_out->t_P, &tAP);

  VSSetMATHeaderWithType(&tAPAT, 4, 4, 1, sizeof(float), (void *)fAPAT, T_FLOAT);
  VSMultMat(&tAP, &pt_out->t_AT, &tAPAT);

  VSAddMat(&tAPAT, &pt_out->t_QT, &pt_out->t_PP);

  return VS_SUCCESS;
}

int RunUpdateKalmanFilter_2D(ALG_KF_2D *pt_out, VSPoint2f *pt_obs) {
  VSMat     tPPHT, t_HXP, tKH, tIKH;
  VSMat     tHPP, t_Xt;
  VSMat     tHPPHT, tHPPHTR, tInv;
  VSMat     tK, tInnov, tInnovK;
  float     fPPHT[4 * 2], fHPP[2 * 4], fHPPHT[2 * 2], fHPPHTR[2 * 2], fInv[2 * 2];
  float     fK[4 * 2], fHXP[2], f_Xt[2], fInnov[2], fInnovK[4], fKH[4 * 4], fIKH[4 * 4];

  /* 1 kalman gain: K = PP*H'*inv(H*PP*H'+R); */
  VSSetMATHeaderWithType(&tPPHT, 2, 4, 1, sizeof(float), (void *)fPPHT, T_FLOAT);
  VSMultMat(&pt_out->t_PP, &pt_out->t_HT, &tPPHT);

  VSSetMATHeaderWithType(&tHPP, 4, 2, 1, sizeof(float), (void *)fHPP, T_FLOAT);
  VSMultMat(&pt_out->t_H, &pt_out->t_PP, &tHPP);

  VSSetMATHeaderWithType(&tHPPHT, 2, 2, 1, sizeof(float), (void *)fHPPHT, T_FLOAT);
  VSMultMat(&tHPP, &pt_out->t_HT, &tHPPHT);

  VSSetMATHeaderWithType(&tHPPHTR, 2, 2, 1, sizeof(float), (void *)fHPPHTR, T_FLOAT);
  VSAddMat(&tHPPHT, &pt_out->t_R, &tHPPHTR);

  VSSetMATHeaderWithType(&tInv, 2, 2, 1, sizeof(float), (void *)fInv, T_FLOAT);
  VSInverseArray(fHPPHTR, fInv, 2);
  //VSInverseMat2x2(&tHPPHTR, &tInv);

  VSSetMATHeaderWithType(&tK, 2, 4, 1, sizeof(float), (void *)fK, T_FLOAT);
  VSMultMat(&tPPHT, &tInv, &tK);

  /* 2 innovative: [cc(t),cr(t)]' - H*xp) */
  VSSetMATHeaderWithType(&t_HXP, 1, 2, 1, sizeof(float), (void *)fHXP, T_FLOAT);
  VSMultMat(&pt_out->t_H, &pt_out->t_XP, &t_HXP);

  VSSetMATHeaderWithType(&t_Xt, 1, 2, 1, sizeof(float), (void *)f_Xt, T_FLOAT);
  f_Xt[0] = pt_obs->x;
  f_Xt[1] = pt_obs->y;

  VSSetMATHeaderWithType(&tInnov, 1, 2, 1, sizeof(float), (void *)fInnov, T_FLOAT);
  VSSubMat(&t_Xt, &t_HXP, &tInnov);

  /* Kalman gain multipled by innovative  */
  VSSetMATHeaderWithType(&tInnovK, 1, 4, 1, sizeof(float), (void *)fInnovK, T_FLOAT);
  VSMultMat(&tK, &tInnov, &tInnovK);

  /* update state vector: Xt = xp + K*([cc(t),cr(t)]' - H*xp)*/
  VSAddMat(&tInnovK, &pt_out->t_XP, &pt_out->t_X);

  /*3. P = (eye(4)-K*H)*PP */
  VSSetMATHeaderWithType(&tKH, 4, 4, 1, sizeof(float), (void *)fKH, T_FLOAT);
  VSMultMat(&tK, &pt_out->t_H, &tKH);

  VSSetMATHeaderWithType(&tIKH, 4, 4, 1, sizeof(float), (void *)fIKH, T_FLOAT);
  VSSubMat(&pt_out->t_I, &tKH, &tIKH);

  VSMultMat(&tIKH, &pt_out->t_PP, &pt_out->t_P);

  return VS_SUCCESS;
}

int BuildKalmanFilter_3D(ALG_KF_3D *pt_out) {
  VSSetMATHeaderWithType(&pt_out->t_X, 1, 3, 1, sizeof(float), (void *)&pt_out->af_X, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_XP, 1, 3, 1, sizeof(float), (void *)&pt_out->af_XP, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_A, 3, 3, 1, sizeof(float), (void *)&pt_out->af_A, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_AT, 3, 3, 1, sizeof(float), (void *)&pt_out->af_AT, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_b, 1, 3, 1, sizeof(float), (void *)&pt_out->af_b, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_P, 3, 3, 1, sizeof(float), (void *)&pt_out->af_P, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_PP, 3, 3, 1, sizeof(float), (void *)&pt_out->af_PP, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_H, 3, 3, 1, sizeof(float), (void *)&pt_out->af_H, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_HT, 3, 3, 1, sizeof(float), (void *)&pt_out->af_HT, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_R, 3, 3, 1, sizeof(float), (void *)&pt_out->af_R, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_Q, 3, 3, 1, sizeof(float), (void *)&pt_out->af_Q, T_FLOAT);
  VSSetMATHeaderWithType(&pt_out->t_QT, 3, 3, 1, sizeof(float), (void *)&pt_out->af_QT, T_FLOAT);

  pt_out->af_I[0] = 1.f; pt_out->af_I[4] = 1.f; pt_out->af_I[8] = 1.f;
  VSSetMATHeaderWithType(&pt_out->t_I, 3, 3, 1, sizeof(float), (void *)&pt_out->af_I, T_FLOAT);

  return VS_SUCCESS;
}

int InitKalmanFilter_3D(ALG_KF_3D *pt_out, VSPoint3f *pt_init_state, float *af_P, float *af_PP) {
  VSPoint3f t_pred, t_cur;

  /* use current observation as prediction */
  pt_out->af_X[0] = pt_init_state->x;
  pt_out->af_X[1] = pt_init_state->y;
  pt_out->af_X[2] = pt_init_state->z;

  pt_out->af_XP[0] = pt_init_state->x;
  pt_out->af_XP[1] = pt_init_state->y;
  pt_out->af_XP[2] = pt_init_state->z;

  memcpy(pt_out->af_P, af_P, sizeof(pt_out->af_P));
  memcpy(pt_out->af_PP, af_PP, sizeof(pt_out->af_PP));

  RunKalmanFilter_3D(pt_out, pt_init_state, &t_pred, &t_cur);

  return VS_SUCCESS;
}

int RunKalmanFilter_3D(ALG_KF_3D *pt_out, VSPoint3f *pt_obs, VSPoint3f *pt_pred, VSPoint3f *pt_cur) {
  /*
  * for convenience for the whole detect(based on prediction)-match-track(based on current detection) framework,
  * at first,kalman tracking is updated for current frame based on prediction at previous moment and observation at
  * current moment,meanwhile currently updated state as outputing model.
  * then, prediction is performed for usage at last frame
  */
  RunUpdateKalmanFilter_3D(pt_out, pt_obs);

  pt_cur->x = pt_out->af_X[0];
  pt_cur->y = pt_out->af_X[1];
  pt_cur->z = pt_out->af_X[2];

  RunPredictionKalmanFilter_3D(pt_out, pt_pred);

  return VS_SUCCESS;
}

int RunPredictionKalmanFilter_3D(ALG_KF_3D *pt_out, VSPoint3f *pt_pred) {
  VSMat         tAX, tAP, tAPAT;
  float         fAX[3], fAP[3 * 3], fAPAT[3 * 3];
  /*
  *  so the predictions include:
  *    Y(n|n-1) = A*Y(n-1|n-1) + b
  *    X(n|n-1) = H*Y(n|n-1)
  *    PP(n|n-1) = A*Pyy()*A' + Q';
  */
  /* the predicted states vector: Y(n|n-1) = A*Y(n-1|n-1) + b */
  VSSetMATHeaderWithType(&tAX, 1, 3, 1, sizeof(float), (void *)fAX, T_FLOAT);
  VSMultMat(&pt_out->t_A, &pt_out->t_X, &tAX);

  VSAddMat(&tAX, &pt_out->t_b, &pt_out->t_XP);
  pt_pred->x = pt_out->af_XP[0];
  pt_pred->y = pt_out->af_XP[1];
  pt_pred->z = pt_out->af_XP[2];

  /* Ry(n|n-1) */
  VSSetMATHeaderWithType(&tAP, 3, 3, 1, sizeof(float), (void *)fAP, T_FLOAT);
  VSMultMat(&pt_out->t_A, &pt_out->t_P, &tAP);

  VSSetMATHeaderWithType(&tAPAT, 3, 3, 1, sizeof(float), (void *)fAPAT, T_FLOAT);
  VSMultMat(&tAP, &pt_out->t_AT, &tAPAT);

  VSAddMat(&tAPAT, &pt_out->t_QT, &pt_out->t_PP);

  return VS_SUCCESS;
}

int RunUpdateKalmanFilter_3D(ALG_KF_3D *pt_out, VSPoint3f *pt_obs) {
  VSMat     tPPHT, t_HXP, tKH, tIKH;
  VSMat     tHPP, t_Xt;
  VSMat     tHPPHT, tHPPHTR, tInv;
  VSMat     tK, tInnov, tInnovK;
  float     fPPHT[3 * 3], fHPP[3 * 3], fHPPHT[3 * 3], fHPPHTR[3 * 3], fInv[3 * 3];
  float     fK[3 * 3], fHXP[3], f_Xt[3], fInnov[3], fInnovK[3], fKH[3 * 3], fIKH[3 * 3];

  /* 1 kalman gain: K = PP*H'*inv(H*PP*H'+R); */
  VSSetMATHeaderWithType(&tPPHT, 3, 3, 1, sizeof(float), (void *)fPPHT, T_FLOAT);
  VSMultMat(&pt_out->t_PP, &pt_out->t_HT, &tPPHT);

  VSSetMATHeaderWithType(&tHPP, 3, 3, 1, sizeof(float), (void *)fHPP, T_FLOAT);
  VSMultMat(&pt_out->t_H, &pt_out->t_PP, &tHPP);

  VSSetMATHeaderWithType(&tHPPHT, 3, 3, 1, sizeof(float), (void *)fHPPHT, T_FLOAT);
  VSMultMat(&tHPP, &pt_out->t_HT, &tHPPHT);

  VSSetMATHeaderWithType(&tHPPHTR, 3, 3, 1, sizeof(float), (void *)fHPPHTR, T_FLOAT);
  VSAddMat(&tHPPHT, &pt_out->t_R, &tHPPHTR);

  VSSetMATHeaderWithType(&tInv, 3, 3, 1, sizeof(float), (void *)fInv, T_FLOAT);
  VSInverseArray(fHPPHTR, fInv, 3);
  //VSInverseMat3x3(&tHPPHTR, &tInv);

  VSSetMATHeaderWithType(&tK, 3, 3, 1, sizeof(float), (void *)fK, T_FLOAT);
  VSMultMat(&tPPHT, &tInv, &tK);

  /* 2 innovative: [cc(t),cr(t)]' - H*xp) */
  VSSetMATHeaderWithType(&t_HXP, 1, 3, 1, sizeof(float), (void *)fHXP, T_FLOAT);
  VSMultMat(&pt_out->t_H, &pt_out->t_XP, &t_HXP);

  VSSetMATHeaderWithType(&t_Xt, 1, 3, 1, sizeof(float), (void *)f_Xt, T_FLOAT);
  f_Xt[0] = pt_obs->x;
  f_Xt[1] = pt_obs->y;
  f_Xt[2] = pt_obs->z;

  VSSetMATHeaderWithType(&tInnov, 1, 3, 1, sizeof(float), (void *)fInnov, T_FLOAT);
  VSSubMat(&t_Xt, &t_HXP, &tInnov);

  /* Kalman gain multipled by innovative  */
  VSSetMATHeaderWithType(&tInnovK, 1, 3, 1, sizeof(float), (void *)fInnovK, T_FLOAT);
  VSMultMat(&tK, &tInnov, &tInnovK);

  /* update state vector: Xt = xp + K*([cc(t),cr(t)]' - H*xp)*/
  VSAddMat(&tInnovK, &pt_out->t_XP, &pt_out->t_X);

  /*3. P = (eye(4)-K*H)*PP */
  VSSetMATHeaderWithType(&tKH, 3, 3, 1, sizeof(float), (void *)fKH, T_FLOAT);
  VSMultMat(&tK, &pt_out->t_H, &tKH);

  VSSetMATHeaderWithType(&tIKH, 3, 3, 1, sizeof(float), (void *)fIKH, T_FLOAT);
  VSSubMat(&pt_out->t_I, &tKH, &tIKH);

  VSMultMat(&tIKH, &pt_out->t_PP, &pt_out->t_P);

  return VS_SUCCESS;
}
