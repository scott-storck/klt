# Test KLT
# S.P.Storck 2017-08-10

import sys
import time
from XMinter import xm, pipe_on, pipe_off
from primitive import m_pick, m_get_switch

def main():
    sim_dur_s = float(m_pick(1, 1.0))
    environ = int(m_pick(2, 2))
    rx_bw_hz = float(m_pick(3, 5e6))
    plots = m_get_switch("PLOTS", 0)

    sim_fname = "jetsim"
    klen = 32
    olap = 0.0
    num_eig = 1

    if 0:
        tstart = time.time()
        xm.jetsim(sim_fname, sim_dur_s, environ, rx_bw_hz)
        print "JETsim took %.3f s"%(time.time() - tstart)
        xm.verify("on")
        tstart = time.time()
        xm.klt(sim_fname, "kltb", None, None, klen, olap, klen, num_eig, olap)
        #dbg="valgrind --tool=callgrind")
        print "KLT took %.3f s"%(time.time() - tstart)
        if num_eig > 1:
            xm.fsvi("kltb(fs=%d)"%(klen), "klt", 0, num_eig, "in")
            xm.headermod("klt", "t=1000")
            xm.detect("klt", None, "fmk")
        else:
            xm.noop("kltb(fs=0)", "klt")
            xm.detect("klt", None, "fmk")
        # Display
        xm.xplot("%s|klt"%(sim_fname), xn="klt %d"%(klen), bg=True, bs=True, auto=3,
                 xl=0, xw=800, xt=0, xh=285)
        xm.detect(sim_fname, None, "fm0")
        xm.xplot("fm0|fmk", xn="klt fmd %d"%(klen), bg=True, bs=True, auto=3,
                 xl=0, xw=800, xt=340, xh=285)
        xm.verify("off")
    elif 1:
        tstart = time.time()
        xm.jetsim(sim_fname, sim_dur_s, environ, rx_bw_hz)
        print "JETsim took %.3f s"%(time.time() - tstart)
        xm.verify("on")
        tstart = time.time()
        xm.kltdet(sim_fname, "fmk", klen, num_eig, "100e9|-50e9")
        print "KLTDET took %.3f s"%(time.time() - tstart)
        xm.thin("fmk", "fmk1", 1, None, 2)
        xm.thin("fmk", "fmk2", 2, None, 2)
        # Display
        xm.xplot("fmk1|fmk2", xn="kltdet fmd %d"%(klen), bg=True, bs=True, auto=3,
                 xl=0, xw=800, xt=00, xh=400)
        xm.verify("off")
    else:
        pipe_on()
        if sim_dur_s <= 0.0:
            xm.xpipemonitor(ncol=8, xl=0, xt=0, xw=480, xh=200)
        pipe_off()
    return 0

if __name__ == "__main__":
    sys.exit(main())
