
int KppDecomp( KPP_REAL *JVS )
{
KPP_REAL W[KPP_NVAR];
KPP_REAL a;
int k, kk, j, jj;

  for( k = 0; k < KPP_NVAR; k++ ) {
    if( JVS[ LU_DIAG[k] ] == 0.0 ) return k+1;
    for( kk = LU_CROW[k]; kk < LU_CROW[k+1]; kk++ )
      W[ LU_ICOL[kk] ] = JVS[kk];
    for( kk = LU_CROW[k]; kk < LU_DIAG[k]; kk++ ) {
      j = LU_ICOL[kk];
      a = -W[j] / JVS[ LU_DIAG[j] ];
      W[j] = -a;
      for( jj = LU_DIAG[j]+1; jj < LU_CROW[j+1]; jj++ )
        W[ LU_ICOL[jj] ] += a*JVS[jj];
    }
    for( kk = LU_CROW[k]; kk < LU_CROW[k+1]; kk++ )
      JVS[kk] = W[ LU_ICOL[kk] ];
  }
  return 0;
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Sparse LU factorization, complex
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int KppDecompCmplxR( KPP_REAL JVSR[], KPP_REAL JVSI[] )
{
   KPP_REAL WR[NVAR], WI[NVAR];
   KPP_REAL ar, ai, den;
   int k, kk, j, jj;

   for( k = 0; k < NVAR; k++ ) {
	if( JVSR[ LU_DIAG[k] ] == 0.0 ) return k+1;
	if( JVSI[ LU_DIAG[k] ] == 0.0 ) return k+1;
	for( kk = LU_CROW[k]; kk < LU_CROW[k+1]; kk++ ) {
		WR[ LU_ICOL[kk] ] = JVSR[kk];
		WI[ LU_ICOL[kk] ] = JVSI[kk];
	}
	for( kk = LU_CROW[k]; kk < LU_DIAG[k]; kk++ ) {
	   j = LU_ICOL[kk];
	   den = JVSR[LU_DIAG[j]]*JVSR[LU_DIAG[j]] + JVSI[LU_DIAG[j]]*JVSI[LU_DIAG[j]];
	   ar = -(WR[j]*JVSR[LU_DIAG[j]] + WI[j]*JVSI[LU_DIAG[j]])/den;
	   ai = -(WI[j]*JVSR[LU_DIAG[j]] - WR[j]*JVSI[LU_DIAG[j]])/den;
	   WR[j] = -ar;
	   WI[j] = -ai;
	   for( jj = LU_DIAG[j]+1; jj < LU_CROW[j+1]; jj++ ) {
		   WR[ LU_ICOL[jj] ] = WR[ LU_ICOL[jj] ] + ar*JVSR[jj] - ai*JVSI[jj];
		   WI[ LU_ICOL[jj] ] = WI[ LU_ICOL[jj] ] + ar*JVSI[jj] + ai*JVSR[jj];
	   }
	}
	for( kk = LU_CROW[k]; kk < LU_CROW[k+1]; kk++ ) {
	   JVSR[kk] = WR[ LU_ICOL[kk] ];
	   JVSI[kk] = WI[ LU_ICOL[kk] ];
	}
   }
   return 0;
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	Complex sparse solve subroutine using indirect addressing
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
void KppSolveCmplxR(KPP_REAL JVSR[], KPP_REAL JVSI[], KPP_REAL XR[], KPP_REAL XI[])
{
   int i, j;
   KPP_REAL sumr, sumi, den;

   for ( i = 0; i < NVAR; i++ ) {
	for ( j = LU_CROW[i]; j < LU_DIAG[i]; j++ ) {
	   XR[i] = XR[i] - (JVSR[j]*XR[LU_ICOL[j]] - JVSI[j]*XI[LU_ICOL[j]]);
	   XI[i] = XI[i] - (JVSR[j]*XI[LU_ICOL[j]] + JVSI[j]*XR[LU_ICOL[j]]);
	}
   }
   
   for ( i = NVAR-1; i >= 0; i-- ) {
	sumr = XR[i];
	sumi = XI[i];
	for ( j = LU_DIAG[i]+1; j < LU_CROW[i+1]; j++) {
	   sumr = sumr - (JVSR[j]*XR[LU_ICOL[j]] - JVSI[j]*XI[LU_ICOL[j]]);
	   sumi = sumi - (JVSR[j]*XI[LU_ICOL[j]] + JVSI[j]*XR[LU_ICOL[j]]);
	}
	den = JVSR[LU_DIAG[i]]*JVSR[LU_DIAG[i]] + JVSI[LU_DIAG[i]]*JVSI[LU_DIAG[i]];
  	XR[i] = (sumr*JVSR[LU_DIAG[i]] + sumi*JVSI[LU_DIAG[i]])/den; 
  	XI[i] = (sumi*JVSR[LU_DIAG[i]] - sumr*JVSI[LU_DIAG[i]])/den;
   }
}
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	END FUNCTION KppSolveCmplxR
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
