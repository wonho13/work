###########################################################################
# program: signal_utilities.py
# author: Tom Irvine
# version: 1.2
# date: November 12, 2013
# 
###########################################################################

from __future__ import print_function


from scipy.signal import lfilter

from numpy import zeros,floor,array
from numpy import sqrt,pi,log

from math import cos,sin,tan

from tompy import read_two_columns_from_dialog

###############################################################################

class MEAN_FILTER:
    
# b=amplitude
# num=number of points    
# np=number of passes
# w=window size
    
    def __init__(self,b,num,np,w):
    
        self.b=b
        self.num=num
        self.np=np
        self.w=w    
        self.np=np
        self.mf=[]
       
    def mean_filter(self):
        k=int(floor(float(self.w-1)/2.))
        last=self.num
        
        trend=zeros(self.num,'f')
        
        for m in range (0,self.np):

            for i in range (0,last):

                ave=0.
                n=0

                for j in range ((i-k),(i+k+1)):  

                    if(j>=0 and j<last ):

                        ave=ave+self.b[j]
                        n=n+1
          	   
                if(n>1):
                    trend[i]=ave/float(n)

        self.mf=self.b-trend
        
        return self.mf

###############################################################################

class BUTTERWORTH:

    def __init__(self,l,f,fh,fl,dt,iband,iphase,y):
        """       
        l=order  (typically 6)
        
        iband
          1= lowpass
          2= highpass 
          3= bandpass
          
        iphase  
          1=yes (appropriate for transient signal)
          2=no  (appropriate for steady-state signal) 
          
        """    
        self.l=l
        self.f=f
        self.freq=f
        self.fh=fh
        self.fl=fl
        self.dt=dt
        self.iband=iband
        self.iphase=iphase
        self.y=y
        
        self.om=0
        
        self.a=zeros((4,4),'f')	
        self.b=zeros((4,4),'f')
        
        self.alpha=zeros(2*self.l,'f')
        
        self.s=(1+1j)*zeros(20,'f')
        
        self.ns=len(y)
        
        self.ik=0
        
        self.yt=zeros(self.ns,'f')        
        
            
    def Butterworth_filter_main(self):
        

        if(self.iband !=3):
            BUTTERWORTH.coefficients(self)
    
        if(self.iband == 1 or self.iband ==2):
            BUTTERWORTH.applymethod(self)

        if(self.iband == 3):
            self.f=self.fh
            self.freq=self.f
            
            #print("\n Step 1")
            self.iband=2
    
            BUTTERWORTH.coefficients(self)
            BUTTERWORTH.applymethod(self)

            self.f=self.fl
            self.freq=self.f

            #print("\n Step 2")
            self.iband=1
    
            BUTTERWORTH.coefficients(self)
            BUTTERWORTH.applymethod(self)   
            
        return self.y        
            
    @classmethod    
    def coefficients(cls,self):
    
        self.a=zeros((4,4),'f')	
        self.b=zeros((4,4),'f')		
    
#*** normalize the frequency ***

        targ=pi*self.f*self.dt   # radians
    
        # print (" targ = %8.4g " %targ)
             
        self.om=tan(targ)   
    
        # print ("   om = %8.4g " %self.om)

#*** solve for the poles *******

        BUTTERWORTH.poles(self)

#*** solve for alpha values ****

        # print("\n alpha ")    
    
        self.alpha=zeros(2*self.l,'f')
        self.alpha=2*self.s.real
    
##    for i in range(0,len(alpha)):
##        print ("  %5.3f +j %5.3f " %(alpha[i].real,alpha[i].imag))

#*** solve for filter coefficients **

        if( self.iband == 1 ):
            BUTTERWORTH.lco(self)
        else:
            BUTTERWORTH.hco(self)
    
#*** plot digital transfer function **

#    dtrans();

#*** check stability ****************
    
        BUTTERWORTH.stab(self)
    
    
    
    @classmethod
    def applymethod(cls,self):
        
        if(self.iphase==1):
            self.apply(self)
            self.apply(self) 
        else:	
            self.apply(self)
        
    

    @classmethod
    def stage1(cls,self):
             
        self.yt=zeros(self.ns,'f')

        bc=self.b[self.ik][0:3]
        ac=self.a[self.ik][0:3]
        ac[0]=1
    
        self.yt=lfilter(bc, ac, self.y, axis=-1, zi=None)      


 
    @classmethod
    def stage2(cls,self):
    
        self.y=zeros(self.ns,'f')
    
        bc=self.b[self.ik][0:3]
        ac=self.a[self.ik][0:3]
        ac[0]=1

        self.y=lfilter(bc, ac, self.yt, axis=-1, zi=None)  
    
        
    @classmethod
    def apply(cls,self):
        
        BUTTERWORTH.coefficients(self)
 
    
        if(self.iphase==1):	

            yr=zeros(self.ns,'f')
            for i in range(0,int(self.ns)):
                yr[self.ns-1-i]=self.y[i]

            self.y=yr
            

#  cascade stage 1

        # print("\n  stage 1")
        self.ik=1
        BUTTERWORTH.stage1(self)

#  cascade stage 2

        # print("  stage 2");
        self.ik=2
        BUTTERWORTH.stage2(self);
       
#  cascade stage 3

        # print("  stage 3");
        self.ik=3
        BUTTERWORTH.stage1(self);
    
        self.y=self.yt

    

    @classmethod	
    def stab(cls,self):
    
        a1=0
        d1=0 
        d2=0 
        d3=0
        dlit=0

        at1=0
        at2=0
        als=0.5e-06
        h2=0

        als*=6.
    
        # print ("\n stability reference threshold= %14.7e " %als)

        for i in range(1,int((self.l/2)+1)):
        
            at1= -self.a[i][1]
            at2= -self.a[i][2]

#       print("\n\n stability coordinates: (%12.7g, %14.7g) ",at1,at2);
        
            h2=at2
 
            a1=h2-1.
            d3=at1-a1
         
            a1=1.-h2
            d2=a1-at1
            d1=at2+1.
		
#       print("\n d1=%14.5g  d2=%14.5g  d3=%14.5g",d1,d2,d3);

            dlit=d1

            if(dlit > d2):
                dlit=d2
            if(dlit > d3):
                dlit=d3

            # print ("\n stage %ld     dlit= %14.5g " %(i, dlit))

            # if(dlit > als):
                # print (" good stability")  			
				
            # if( (dlit < als) and (dlit > 0.)):		  
                # print(" marginally unstable ");
            
            # if(dlit < 0.):
                # print (" unstable ")	  	
                # print ("\n")

################################################################################  
    @classmethod	
    def lco(cls,self):
    
        om2=self.om**2

        for k in range(1,int((self.l/2)+1)):
    
            den = om2-self.alpha[k-1]*self.om+1.
		
            self.a[k][0]=0.
            self.a[k][1]=2.*(om2 -1.)/den
            self.a[k][2]=( om2 +self.alpha[k-1]*self.om+ 1.)/den

            self.b[k][0]=om2/den
            self.b[k][1]=2.*self.b[k][0]
            self.b[k][2]=self.b[k][0]

            # print ("\n filter coefficients")		
            # print (" a[%i][1]=%10.5g  a[%i][2]=%10.5g" \
            #                                   %(k,self.a[k][1],k,self.a[k][2]))
            # print (" b[%i][0]=%10.5g  b[%i][1]=%10.5g  b[%i][2]=%10.5g" \
            #                    %(k,self.b[k][0],k,self.b[k][1],k,self.b[k][2]))
    
        # print ("\n")
        

################################################################################  
    @classmethod	
    def hco(cls,self):
    
        # print ("\n filter coefficients")
    
        om2=self.om**2

        for k in range(1,int((self.l/2)+1)):
    
            den = om2-self.alpha[k-1]*self.om+1.    
		
            self.a[k][0]=0.
            self.a[k][1]=2.*(-1.+ om2)/den
            self.a[k][2]=( 1.+self.alpha[k-1]*self.om+ om2)/den

            self.b[k][0]= 1./den;
            self.b[k][1]=-2.*self.b[k][0]
            self.b[k][2]=    self.b[k][0]
        
            # print ("\n a[%i][1]=%10.5g  a[%i][2]=%10.5g" \
            #                                  %(k,self.a[k][1],k,self.a[k][2]))
            # print (" b[%i][0]=%10.5g  b[%i][1]=%10.5g  b[%i][2]=%10.5g" \
            #                   %(k,self.b[k][0],k,self.b[k][1],k,self.b[k][2]))
            # print ("\n")
        

################################################################################  
    @classmethod	
    def poles(cls,self):
        arg=0
        a1=0
        a2=complex(0.,0.)
        h=complex(0.,0.)
        theta=complex(0.,0.)
    
        self.s=(1+1j)*zeros(20,'f')    
    
#    print("\n  calculate print ");

        # print ("\n poles ")
	
        for k in range(0,int(2*self.l)):
            arg=(2.*(k+1) +self.l-1)*pi/(2.*self.l)
            self.s[k]=cos(arg)+sin(arg)*(1j)
            # print (" %4.3f  +j %4.3f " %(self.s[k].real,self.s[k].imag))
         
        for i in range(0,201):   
            arg = i/40.        
        
            h=complex( self.s[0].real,( arg - self.s[0].imag  ))

        for j in range(1,int(self.l)):
            
            theta=complex( -self.s[j].real,( arg - self.s[j].imag ))
            
            temp=h*theta
            h=temp
               
            x=1/h
            h=x
               
            a1 = self.freq*arg
	   
            a2=abs(h)            
           
            a3 = a2**2
         
###############################################################################

def EnterPSD():
    """
    a = frequency column
    b = PSD column
    num = number of coordinates
    slope = slope between coordinate pairs    
    """
    print (" ")
    print (" The input file must have two columns: \n freq(Hz) & accel(G^2/Hz)")
    print (" (Find dialog box) ")

    a,b,num =read_two_columns_from_dialog('Select PSD file')

    print ("\n samples = %d " % num)

    a=array(a)
    b=array(b)
    

    nm1=num-1

    slope =zeros(nm1,'f')


    ra=0

    for i in range (0,nm1):
#
        s=log(b[i+1]/b[i])/log(a[i+1]/a[i])
        
        slope[i]=s
#
        if s < -1.0001 or s > -0.9999:
            ra+= ( b[i+1] * a[i+1]- b[i]*a[i])/( s+1.)
        else:
            ra+= b[i]*a[i]*log( a[i+1]/a[i])

    omega=2*pi*a

    bv=zeros(num,'f') 
    bd=zeros(num,'f') 
        
    for i in range (0,num):         
        bv[i]=b[i]/omega[i]**2
     
    bv=bv*386**2
    rv=0

    for i in range (0,nm1):
#
        s=log(bv[i+1]/bv[i])/log(a[i+1]/a[i])
#
        if s < -1.0001 or s > -0.9999:
            rv+= ( bv[i+1] * a[i+1]- bv[i]*a[i])/( s+1.)
        else:
            rv+= bv[i]*a[i]*log( a[i+1]/a[i])         
         
        
    for i in range (0,num):         
        bd[i]=bv[i]/omega[i]**2
     
    rd=0

    for i in range (0,nm1):
#
        s=log(bd[i+1]/bd[i])/log(a[i+1]/a[i])
#
        if s < -1.0001 or s > -0.9999:
            rd+= ( bd[i+1] * a[i+1]- bd[i]*a[i])/( s+1.)
        else:
            rd+= bd[i]*a[i]*log( a[i+1]/a[i])         


    rms=sqrt(ra)
    three_rms=3*rms
    
    print (" ")
    print (" *** Input PSD *** ")
    print (" ")
 
    print (" Acceleration ")
    print ("   Overall = %10.3g GRMS" % rms)
    print ("           = %10.3g 3-sigma" % three_rms)

    grms=rms

    vrms=sqrt(rv)
    vthree_rms=3*vrms

    print (" ")
    print (" Velocity ") 
    print ("   Overall = %10.3g in/sec rms" % vrms)
    print ("           = %10.3g in/sec 3-sigma" % vthree_rms)

    drms=sqrt(rd)
    dthree_rms=3*drms

    print (" ")
    print (" Displacement ") 
    print ("   Overall = %10.3g in rms" % drms)
    print ("           = %10.3g in 3-sigma" % dthree_rms)

    return a,b,grms,num,slope         