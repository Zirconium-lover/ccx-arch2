/*     CalculiX - damage/fracture extension                              */
/*     dammat.c: what kind of damage material is this?                   */

/* Why this module exists
   ----------------------
   One question - "does this material carry a progressive damage law, and
   which one" - is asked by code at two different levels of the solver.  The
   erosion rules ask it to decide whether an element is even a candidate for
   deletion; the topology census asks it to decide whether a degradation
   value needs the DM2.0 unit shift before it can be compared.

   The answer lived in erosion.c, because erosion was the caller that got
   written first.  That file already carried the right instinct in a comment
   - the constant next to this function was renamed "under the name of the
   question it answers rather than of the one caller that happened to be
   written first" - but the instinct was applied to the NAME and not to the
   FILE, and the function stayed where its first caller was.

   That left an upward call: topology.c, which describes what the model IS,
   reaching into erosion.c, which decides what leaves it.  A model that has
   to ask a policy a question cannot be understood, tested or changed without
   the policy, and the direction of the dependency says the wrong thing about
   which of the two is the more fundamental.  tools/arch.py reports it as a
   layer violation, and this file is the fix: the classifier moves below both
   callers, where a pure statement about a material card belongs.

   What this is
   ------------
   A classifier of the MATERIAL CARD and nothing else.  It reads ndmcon and
   dmcon - the constants the deck declared - and returns a verdict about
   their shape.  It is a pure function: same arguments, same answer, no
   state, no allocation, no output.

   What this is NOT
   ----------------
   Not a reading of damage.  What D is at an integration point is
   calcdamage's answer; this file never looks at it.  Not a decision: it says
   what a material IS, never what should be done about it.  Both of those
   belong to the callers, and keeping them out is what stops this file
   growing back into the one it came from.                               */

#include <stdio.h>
#include "CalculiX.h"

/* Progressive damage material classifier shared by DE1 and DM2.0.
   Rice-Tracey + Evolution=Displacement keeps the historical four-constant
   signature.  DM2.0 is identified by model type 3 and a variable-length
   constant count 3+2*NPOINTS (NPOINTS>=2). */
ITG damage_progressive_material(ITG imat,const ITG *ndmcon,
                                const double *dmcon,ITG ndmat,
                                ITG ntmat)
{
  ITG nconst,type,off;

  if(imat<1) return 0;
  nconst=ndmcon[2*(imat-1)];
  if((dmcon==NULL)||(ndmat<1)||(ntmat<1)) return 0;

  off=1+(ndmat+1)*ntmat*(imat-1);
  type=(ITG)dmcon[off];
  if((type==1)&&(nconst==4)) return 1;
  if((type==3)&&(nconst>=7)&&(((nconst-3)%2)==0)) return 1;

  return 0;
}

static void dammat_chki(const char *what,ITG got,ITG want,ITG *nbad)
{
  if(got!=want) (*nbad)++;
  printf("[DAMMAT]   %-52s %6" ITGFORMAT " (want %" ITGFORMAT ") %s\n",
         what,got,want,(got==want)?"PASS":"FAIL");
}

/* The two cases that decide whether an element is a deletion candidate at
   all.  They came over from erosion_selftest() with the function; three
   materials rather than two because the second case asks about material 3,
   and both tables are indexed by material number - ndmcon at 2*(imat-1),
   dmcon at 1+(ndmat+1)*ntmat*(imat-1), which is 11 for imat=3.  Sized for
   two, the question read off the end of both arrays; -Wall said so and it
   was right.

   The DM2.0 cases are new here.  The classifier has always accepted type 3
   with an odd-shaped constant count and nothing tested that it rejects the
   shapes next to it, which is the half of a classifier that is easy to get
   wrong: 3+2*NPOINTS with NPOINTS>=2 means seven, nine, eleven - and six and
   eight must come back zero. */
ITG dammat_selftest(void)
{
  ITG nbad=0,i;
  ITG ndmcon[8]={4,0,4,0,0,0,7,0};
  double dmcon[20];

  /* off=1+(ndmat+1)*ntmat*(imat-1) with ndmat=4, ntmat=1, so the type of
     material n sits at dmcon[1+5*(n-1)]: 1, 6, 11, 16.  Writing material
     4's type to 11 instead of 16 is not a test that fails, it is a test
     that cannot fail - type came back 0, no type-3 branch was ever
     reached, and all three DM2.0 cases agreed with each other for the
     wrong reason.  Seen only by running it. */
  for(i=0;i<20;i++) dmcon[i]=0.;
  dmcon[1]=1.;                        /* material 1: Rice-Tracey, type 1 */
  dmcon[6]=1.;                        /* material 2: same                */
  dmcon[16]=3.;                       /* material 4: DM2.0, type 3       */

  dammat_chki("progressive material is recognised",
              damage_progressive_material(1,ndmcon,dmcon,4,1),1,&nbad);
  dammat_chki("a material with no damage record is not",
              damage_progressive_material(3,ndmcon,dmcon,4,1),0,&nbad);
  /* The imat<1 guard is the one thing here a value test cannot hold down.
     Its job is to stop ndmcon[2*(imat-1)] reading backwards off the front
     of the array; removing it does not produce a WRONG ANSWER, it produces
     undefined behaviour that happened to return zero on this compiler.
     Mutating `imat<1' to `imat<0' was seen leaving this test green.  The
     assertion below is still worth having - it pins the contract - but the
     guard itself wants a sanitiser build, not a comparison. */
  dammat_chki("material number below one is not",
              damage_progressive_material(0,ndmcon,dmcon,4,1),0,&nbad);
  dammat_chki("a NULL constant table is not",
              damage_progressive_material(1,ndmcon,NULL,4,1),0,&nbad);
  /* ndmat and ntmat size the stride into dmcon.  Zero for either means
     there is no constant table to index, whatever ndmcon claims; nothing
     asked, and relaxing both guards to <0 went unnoticed. */
  dammat_chki("ndmat of zero is not",
              damage_progressive_material(1,ndmcon,dmcon,0,1),0,&nbad);
  dammat_chki("ntmat of zero is not",
              damage_progressive_material(1,ndmcon,dmcon,4,0),0,&nbad);
  dammat_chki("DM2.0 type 3 with 3+2*NPOINTS constants is recognised",
              damage_progressive_material(4,ndmcon,dmcon,4,1),1,&nbad);
  /* The type-1 branch wants EXACTLY four constants - Rice-Tracey plus
     Evolution=Displacement, the historical signature.  Nothing asked what
     happens to a type-1 card with five, and relaxing the test to nconst>=4
     went unnoticed; this is that case. */
  ndmcon[2]=5;
  dammat_chki("type 1 with more than four constants is not",
              damage_progressive_material(2,ndmcon,dmcon,4,1),0,&nbad);
  ndmcon[2]=4;
  dammat_chki("type 1 with exactly four constants is",
              damage_progressive_material(2,ndmcon,dmcon,4,1),1,&nbad);

  /* NPOINTS=1, i.e. nconst=5.  This is the case that actually tests the
     lower bound: the parity rule already forces nconst odd, so >=7 and >=6
     accept exactly the same integers and a mutation between them is not a
     mutation at all.  Five is the first odd count below the bound, and it
     is the shape a DM2.0 card with a single point would have. */
  ndmcon[6]=5;
  dammat_chki("type 3 with one point (nconst=5) is not",
              damage_progressive_material(4,ndmcon,dmcon,4,1),0,&nbad);
  ndmcon[6]=6;
  dammat_chki("type 3 with too few constants is not",
              damage_progressive_material(4,ndmcon,dmcon,4,1),0,&nbad);
  ndmcon[6]=8;
  dammat_chki("type 3 with an even constant count is not",
              damage_progressive_material(4,ndmcon,dmcon,4,1),0,&nbad);

  printf("[DAMMAT] self test: %" ITGFORMAT " failure(s) -- %s\n",
         nbad,nbad?"FAILED":"PASSED");
  fflush(stdout);
  return nbad;
}
