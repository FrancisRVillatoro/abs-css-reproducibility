import json,sys,numpy as np
from scipy.signal import resample

def load(f):
 d=json.load(open(f)); return np.array(d['u'],float),d

def common(a,b,n=3600): return resample(a,n),resample(b,n)
def align_diff(a,b):
 # shift b to a via periodic FFT phase after integer+parabolic lag, return relative norm
 n=len(a); A=np.fft.fft(a); B=np.fft.fft(b); cc=np.fft.ifft(A*np.conj(B)).real; j=int(np.argmax(cc)); jm=(j-1)%n;jp=(j+1)%n; den=cc[jm]-2*cc[j]+cc[jp]; delta=0 if den==0 else .5*(cc[jm]-cc[jp])/den; delta=float(np.clip(delta,-.5,.5)); lag=j+delta; 
 if lag>n/2: lag-=n
 k=2*np.pi*np.fft.fftfreq(n); bs=np.fft.ifft(B*np.exp(-1j*k*lag)).real
 return float(np.linalg.norm(a-bs)/np.linalg.norm(a)),float(lag)
if __name__=='__main__':
 a,da=load(sys.argv[1]);b,db=load(sys.argv[2]);a,b=common(a,b); print(align_diff(a,b))
