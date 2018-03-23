c
c---- tags used for MPI messages
c
      integer TAG_ERROR,TAG_OK,TAG_INPUT,TAG_DIMS
      integer TAG_IMAP,TAG_PATCH,TAG_GRID,TAG_RSTRT
      integer TAG_DIST,TAG_VOL,TAG_FLOW,TAG_FLOW2,TAG_TURB,TAG_TURB2
      integer TAG_DTMIN,TAG_SAVE,TAG_CONV,TAG_TIME,TAG_GRP
c
      parameter (TAG_ERROR =  0)
      parameter (TAG_OK    =  1)
      parameter (TAG_INPUT =  2)
      parameter (TAG_DIMS  =  3)
      parameter (TAG_IMAP  =  4)
      parameter (TAG_PATCH =  5)
      parameter (TAG_GRID  =  6)
      parameter (TAG_RSTRT =  7)
      parameter (TAG_DIST  =  8)
      parameter (TAG_VOL   =  9)
c TAG_FLOW is 100 + iseg (see bcflow.F)
c TAG_FLOW2 is 2100 + iseg (see bcflow.F)
      parameter (TAG_FLOW  = 100)
      parameter (TAG_FLOW2 = 2100)
      parameter (TAG_TURB  = 11)
      parameter (TAG_TURB2 = 4100)
      parameter (TAG_DTMIN = 12)
      parameter (TAG_SAVE  = 13)
      parameter (TAG_CONV  = 14)
      parameter (TAG_TIME  = 15)
      parameter (TAG_GRP   = 16)
