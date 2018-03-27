# pipet

PiPET is a simple, yet completely functional PET acquisition system. Its performances are below the state of the art, it features only two small detectors, but apart from that, it is not too different from any common commercial PET system. Most of its electronics are exposed, so that students can probe the main signals involved during open the photons acquisition. The only parts that are not exposed to the student are the scintillating crystal, the photomultiplier (PMT) and the first layer of position encoding, which is directly coupled to the PMT. These parts are inserted in an opaque plastic case that shields the detector from environmental light (Figure 1).

<p align="center"><img src="https://github.com/gsportelli/pipet/blob/master/img/lyso-in-case.png?raw=true" alt="LYSO matrix"/>Figure 1: Picture of a LYSO matrix inserted in its plastic case. On top of the matrix there are five optical grease drops that will serve for optical coupling with the photomultiplier.
</p>

The architecture of the PiPET acquisition system is shown in (Figure 2). The system is composed of two detector blocks, a FPGA-based motherboard that serves as a coincidence data acquisition processor and motherboard, an adapter card and two cards of the digitization of Anger coordinates, each of which is controlled by an FPGA. The motherboard is finally connected via USB to a PC (raspberry pi 2) that controls the functions through a web interface.

<p align="center"><img src="https://github.com/gsportelli/pipet/blob/master/img/pipet-arch.png?raw=true" alt="pipet architecture"/>Figure 2: Schematic diagram of the PiPET acquisition system.</p>
