from ROOT import gSystem
import os

wcsimdir = 'WCSIMDIR'
wcsim = os.getenv(wcsimdir)
print(wcsimdir)
print(wcsim)
gSystem.Load(wcsim+'/libWCSimRoot.so')
