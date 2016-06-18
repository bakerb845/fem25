          real*8 sig,del,gam,fmom,sol,tconst
          real*8 xlength
          common/sources/sol(6),sig,del,gam,fmom,tconst,xlength
 
          real*8 d0,th,ph,sec
          integer*4 jy,jd,jh,jm
          common/dhead/d0,th,ph,sec, jy,jd,jh,jm

          real*8 to,dist,t1,p1
          integer*4 ivel,mdmin,mdmx,nsdp,init
          common/vary/to,dist,t1,p1, ivel,mdmin,mdmx,nsdp,init

          real*8 sdep
          integer*4 npts,nsrce,jcom,nbran,IRAY,ILOV,nom
          common/same/sdep(lsd),npts,nsrce,jcom,nbran,IRAY,ILOV,nom

          complex*16 resp
          common/resp/resp(lgrm/2+1)

          character  sta*4,chn*4,typ*4
          real*8 ss,dt,df
          integer*4 nscan,iy,id,ih,im
          common/myhed/ss,dt,df, sta,chn,typ,nscan,iy,id,ih,im

