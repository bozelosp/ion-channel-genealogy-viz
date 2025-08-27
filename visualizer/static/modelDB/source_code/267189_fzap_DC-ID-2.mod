: $Id: fzap_DC.mod,v 1.1 2009/04/21 19:17:15 hines Exp hines $

COMMENT
fzap_DC.mod (based on fzap.mod)

A bogus point process that contains the variable x, 
which oscillates starting at t = del >= 0.
Direct current stimulation,
where both del and dur are > 0.

fzap_DC uses the event delivery system to ensure compatibility with adaptive integration.

=================
NOTES AND CAVEATS
=================

1.  If x were a RANGE variable, an assignment statement would 
have to be inserted into proc advance() in order for the 
value of x to be used by other mechanisms--e.g.
proc advance() {
  is_xtra = Fzap_DC[0].x
  fadvance()
}
However, that would be incompatible with adaptive integration.
To eliminate the need for such an assignment statement, x is a 
POINTER.  This preserves compatibility with adaptive integration.

2.  On every fadvance, the statements that evaluate Fzap_DC's x 
should be executed before the statements in any client mechanism 
that relies on the value of Fzap_DC's x.  To that end, the value of 
x is computed in a BEFORE BREAKPOINT block, which will take care
of any client mechanism that uses Fzap_DC's x in a BREAKPOINT block.

However, some client mechanisms may have their own 
BEFORE BREAKPOINT blocks that need the value of Fzap_DC's x.  
xtra is such a mechanism.  In this situation, care is required 
to ensure that the statements in Fzap_DC's BEFORE BREAKPOINT block
are executed first.  This can be done by compiling the mod file 
that defines Fzap_DC _before_ the client mechanism's mod file.

There are two ways to make this happen:
A.  Invoke nrnivmodl with a command line that presents the file 
names in the desired sequence.  UNIX/Linux users may be quite 
comfortable with this.
B.  Choose mod file names so that Fzap_DC's mod file appears before 
the name of any client mod files in an alphabetical listing.
For the example of Fzap_DC and xtra, the file names fzap_DC.mod and 
xtra.mod would be quite suitable.  This is more convenient for 
users of all operating systems, but especially MSWin and OS X, 
whose users are accustomed to compiling all mod files in a 
directory with mknrndll or "drag and drop," respectively.

12/11/2008 NTC
ENDCOMMENT

NEURON {
  POINT_PROCESS Fzap_DC
  RANGE del, dur, amp
  POINTER x
}

PARAMETER {
  del (ms)
  dur (ms)
  amp (1)
}

ASSIGNED {
  x (1)
  on (1)
}

INITIAL {
  x = 0
  on = 0

  if (del<0) { del=0 }
  if (dur<0) { dur=0 }

  : do nothing if dur == 0
  if (dur>0) {
    net_send(del, 1)  : to turn it on
    net_send(del+dur, 1)  : to stop
  }
}

BEFORE BREAKPOINT {
  if (on==0) {
    x = 0
  } else {
    x = amp
  }
}

NET_RECEIVE (w) {
  : respond only to self-events with flag > 0
  if (flag == 1) {
VERBATIM
	on = (double)(on == 0.0);
ENDVERBATIM
  }
}
