void rstrct(double ***uc, double ***uf, int nc)
/* 
   Half-weighting restriction. nc is the coarse-grid dimension. The fine-grid solution is input in
   uf[1..2*nc-1][1..2*nc-1], the coarse-grid solution is returned in uc[1..nc][1..nc].
*/
{
  int ic,iif,jc,jf,kc,kf,ncc=2*nc-1;
  /* Interior points.*/
  for (kf=3,kc=2;kc<nc;kc++,kf+=2) {
    for (jf=3,jc=2;jc<nc;jc++,jf+=2) { 
      for (iif=3,ic=2;ic<nc;ic++,iif+=2) {
         uc[ic][jc][kc]=0.50*uf[iif][jf][kf]+0.08333333333333*(uf[iif+1][jf][kf]+uf[iif-1][jf][kf]
					+uf[iif][jf+1][kf]+uf[iif][jf-1][kf]
					+uf[iif][jf][kf+1]+uf[iif][jf][kf-1]);
      }
    }
  }
  /* Boundary planes. */
  for (iif=1,ic=1;ic<=nc;ic++,iif+=2) { 
    for (kc=1,kf=1;kc<=nc;kc++,kf+=2) {
       uc[ic][1][kc]=uf[iif][1][kf];
       uc[ic][nc][kc]=uf[iif][ncc][kf];
    }
  }
  for (iif=1,ic=1;ic<=nc;ic++,iif+=2) {
    for (kc=1,kf=1;kc<=nc;kc++,kf+=2) {
       uc[1][ic][kc]=uf[1][iif][kf];
       uc[nc][ic][kc]=uf[ncc][iif][kf];
    }
  }
  for (iif=1,ic=1;ic<=nc;ic++,iif+=2) {
    for (kc=1,kf=1;kc<=nc;kc++,kf+=2) {
       uc[ic][kc][1]=uf[iif][kf][1];
       uc[ic][kc][nc]=uf[iif][kf][ncc];
    }
  }
}
