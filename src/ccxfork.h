/*     CalculiX - damage/fracture extension                              */
/*     ccxfork.h: everything this fork adds to the nonlinear solve.      */

/* Why this file exists
   --------------------
   These 1400 lines used to sit inside CalculiX.h, between the stock
   declarations of splitline() and spooles().  CalculiX.h is included by 210
   of the 211 .c files in this tree, so changing ANY extension interface -
   adding an argument to one function nobody else calls - recompiled the
   whole tree.  That is the price of every experiment, paid by whoever tries
   anything, and it is not a question of style.

   Moving the block here does not make the declarations any better.  The
   widest of them still takes twenty-two arguments.  What it changes is who
   pays: thirty files include this, not two hundred and ten.

   This is deliberately ONE file rather than twenty-three.  The partition by
   owning module is written and correct - tools/arch.py already knows which
   module owns every symbol - but it is a twenty-three-way split with a long
   tail of compile errors, and the rule this work follows is that every step
   leaves the tree building and the gate green.  So the cheap, safe 85% goes
   in first and the split follows module by module, each with its own commit
   and its own A/B.  ccxopt.h and logview.h are the first two out: they are
   the platform, they almost never change, and they are what the twenty-eight
   files that only wanted ccxopt_getenv() were paying two hundred and ten
   files' worth of rebuild for.                                          */

#ifndef CCX_FORK_H
#define CCX_FORK_H

#include "CalculiX.h"

/* The contexts are forward-declared here so that the order of declarations
   in the rest of this file does not matter.  Without it, a function that
   takes a trialctx has to be declared BELOW trialctx, which makes the
   header's layout a constraint on which conversions can be done next -
   exactly the kind of accidental coupling this work is removing. */
typedef struct trialctx trialctx;
typedef struct nlstate nlstate;
typedef struct erosion_batch erosion_batch;
typedef struct erosion_policy erosion_policy;
typedef struct topo_txn topo_txn;

/* dissipation-based path following; see pathfollow.c for the derivation
   and for the numerical verification of both rows of the bordered system.
   The driver holds no copy of the displacement: the caller projects f_hat
   against the model and passes the two scalars in. */

double pathfollow_dg(double Pn,double lamn,double P,double lam);
void pathfollow_dgrad(double lamn,double Pn,double fr,double ff,
                      double *adur,double *aduf,double *bb);
ITG pathfollow_dlam(double g,double lamn,double Pn,double fr,double ff,
                    double dlmax,double *dlam,ITG *reason);
ITG pathfollow_selftest(void);
ITG pathfollow_legacycheck(void);

ITG pathfollow_arm(double tau,ITG neq);
void pathfollow_disarm(void);
ITG pathfollow_armed(void);
ITG pathfollow_resize(ITG neq);
double pathfollow_lamn(void);
double pathfollow_Pn(void);
double pathfollow_ff(void);
double pathfollow_gettau(void);
void pathfollow_settau(double tau);
ITG pathfollow_refusals(void);
ITG pathfollow_have(void);
const double *pathfollow_fhat(void);
void pathfollow_incstart(void);
void pathfollow_freeze(void);
void pathfollow_unfreeze(void);
ITG pathfollow_frozen(void);
void pathfollow_capture(const double *b,double dlampred);
void pathfollow_setPn(double Pn);
void pathfollow_measure(const double *uf);
ITG pathfollow_project_lambda(double pdu,double *lam);
ITG pathfollow_predictor(double *dlam);
ITG pathfollow_step(double *b,const double *uf,double pdu,double *lam,
                    double dlmax,double *dgout,double *gout,double *dlamout,
                    ITG *reason);
void pathfollow_commit(double lam,double pdu,double *dgcommit);

/* crack-opening control: linear functional of u, exact derivatives */
ITG pathfollow_cod_arm(const double *c,ITG neq);
ITG pathfollow_cod(void);
double pathfollow_cod_target(void);
void pathfollow_cod_settarget(double t);
const double *pathfollow_cod_c(void);
ITG pathfollow_cod_step(double *b,const double *uf,double cu,double *lam,
                        double dlmax,double *gout,double *dlamout,
                        ITG *reason);

/* ---- mixed-mode crack control (crackcontrol.c) ---------------------
   The control functional for a UC6 cohesive interface: a frozen
   linearisation of the effective separation
       deff^2 = max(dn,0)^2 + beta*(ds1^2+ds2^2)
   so that the constraint stays affine during one corrector.  See the
   block comment in crackcontrol.c. */

typedef struct{
  ITG nfacet;      /* UC6 elements in the model                        */
  ITG ndead;       /* of those, deleted (ipkon<0)                      */
  ITG nip;         /* live integration points visited                  */
  ITG ninit;       /* dmax > d0  (initiated)                           */
  ITG nzone;       /* d0 < dmax < df  (process zone)                   */
  ITG nload;       /* of those, currently loading (deff >= dmax)       */
  ITG nfail;       /* dmax >= df (fully failed)                        */
  double area;     /* total live facet area                            */
  double zonearea; /* facet area carrying the process zone             */
  double weight;   /* total weight the functional carries              */
  double shearfrac;/* area-weighted beta*|ds|^2/deff^2 in the zone     */
  double deffmax;  /* largest effective separation                     */
  double dmaxmax;  /* largest committed maximum separation             */
}crackcontrol_census;



/* ---- does the assembled tangent differ the residual? (opcheck.c) ----

   The check for the named hole in the diagnostics, which turned
   out to be already in the tree behind three undeclared switches.  What is
   new here is the DISCRIMINATOR: a central difference converges to the mean
   of the one-sided derivatives at a kink, so its h-plateau cannot tell a
   kink from a wrong tangent.  The two one-sided differences can. */

#define OPCHECK_OK    0   /* one-sided differences agree, and with K_asm   */
#define OPCHECK_KINK  1   /* they disagree, and K_asm is on one branch     */
#define OPCHECK_WRONG 2   /* they agree, and K_asm is on neither           */
#define OPCHECK_BOTH  3   /* they disagree and K_asm is on neither branch  */

typedef struct{
  ITG n;             /* coefficients compared in this column             */
  ITG nok,nkink,nwrong,nboth;
  double scale;      /* largest |K| in the column: everything is relative */
  double relctr;     /* max |K_ctr - K_asm| / scale                       */
  double relside;    /* max |K_fwd - K_bwd| / scale - the kink measure    */
  double relbest;    /* max over rows of min(|K_fwd-K_asm|,|K_bwd-K_asm|) */
  ITG worst,worstdir;/* where relctr peaked                               */
  double wfwd,wbwd,wasm;
  double tol;
}opcheck;

ITG  opcheck_classify(double kfwd,double kbwd,double kasm,double scale,
                      double tol);
void opcheck_begin(opcheck *o,double tol);
void opcheck_column(opcheck *o,ITG nk,ITG mt,const ITG *nactdof,
                    const double *f0,const double *fp,const double *fm,
                    double h,ITG col,
                    const double *ad,const double *au,
                    const ITG *jq,const ITG *irow,const ITG *nzs,ITG nasym);
ITG  opcheck_selftest(void);
void monitor_opcheck(const opcheck *o,ITG iinc,ITG iit,ITG node,ITG dir,
                     double h);
void monitor_opcheck_total(const opcheck *o,ITG iinc,ITG iit,ITG elem,
                           double dam,double h,ITG ncol);

/* ---- the stiffness census (census.c) and its presenter (monitor.c) ----

   How far the assembled nodal stiffness has fallen, as a fraction of each
   node's OWN first positive value.  This loop was written out THREE TIMES,
   character for character, in nonlingeo.c, because it had no owner.  It
   measures; it decides nothing - damstate.c owns the judgement it is read
   against. */

typedef struct{
  ITG below1;        /* nodes below 1e-3 of their own intact diagonal  */
  ITG below2;        /* below 1e-2                                     */
  ITG below3;        /* below 1e-1                                     */
  ITG nonpositive;   /* diagonals that are not positive at all         */
  ITG worst;         /* the extreme node, 1-based, or -1               */
  ITG nnode;         /* nodes that had an intact reference to compare  */
  double worstratio; /* and its ratio                                  */
}stiffcensus;

void stiffcensus_take(stiffcensus *c,ITG nk,const double *addiag,
                      const double *addiag0,const ITG *addok);
ITG  stiffcensus_selftest(void);
void monitor_stiffness(const stiffcensus *c,ITG iinc,double time);



/* ---- THE judgement about what no longer carries load (damstate.c) -----

   One owner for a decision that used to be made in three places for three
   consumers, with the solve never told at all.  See damstate.c. */

typedef struct{
  ITG nk;            /* number of nodes                                   */
  ITG *ok;           /* per node: all three translational dofs are active  */
  ITG *dead;         /* per node: judged no longer load-carrying          */
  double *diag;      /* per node: min assembled diagonal over its dofs     */
  double *diag0;     /* per node: first positive diagonal ever seen        */
  ITG ndead;
  double g;          /* dead below this fraction of the node own intact    */
  ITG allowneg;      /* count a non-positive diagonal as dead              */
}damstate;

void damstate_init(damstate *s,ITG nk,double g,ITG allowneg);
void damstate_free(damstate *s);
void damstate_update(damstate *s,const double *ad,const ITG *nactdof,ITG mt);
ITG  damstate_dead(const damstate *s,ITG node);
ITG  damstate_facet_dead(const double *xstate,ITG nstate,ITG mi0,
                         ITG elem,ITG nip);
ITG  damstate_facet_support(const ITG *ipkon,const char *lakon,const ITG *kon,
                            ITG ne,ITG nk,ITG *nfac);
ITG  damstate_selftest(void);

/* ---- the numbers the convergence judgement is made from (converge.c) --

   Owns the reduction of the solution
   vector into ram/ram1/ram2, the uam high-water mark, the qam running
   average and its floor, and which dofs are allowed to contribute.  It
   decides nothing and prints nothing: checkconvergence.c still owns the
   verdict.  What it EXCLUDED is part of its state, because the promise that
   the excluded residual is reported next to the criterion belongs to the
   thing that does the excluding. */

typedef struct{
  double qam_floor;  /* CCX_DAMAGE_QAM_FLOOR; 0 disables the floor        */
  ITG mask_force;    /* CCX_DAMAGE_AUTOSPC_FORCE armed                    */
  const ITG *mask;   /* per-node load-path judgement (damstate.dead)      */
  ITG mask_nk;
  double qam_peak;   /* largest qam ever seen, for the floor              */
  double excl_max;   /* largest residual excluded by the mask, last call  */
  ITG excl_node;     /* where it was, 1-based                             */
  ITG excl_count;    /* how many dofs were excluded                       */
}converge;

void converge_init(converge *c,double qam_floor,ITG mask_force,
                   const ITG *mask,ITG mask_nk);
void converge_norms(converge *c,const double *b,const ITG *neq,
                    const ITG *nactdofinv,ITG mt,ITG ithermal,ITG mortar,
                    ITG ne,ITG ne0,ITG neold,
                    const double *qa,const double *qamold,ITG jnz,
                    double qau,double ea,
                    double *ram,double *ram1,double *ram2,
                    const double *cam,double *uam,double *qam);
/* Ten arguments became five: the five norms are nlstate's, ithermal is
   trialctx's and mt is mi[1]+1. */
void converge_report(const converge *c,const trialctx *mdl,const nlstate *n,
                     const ITG *nactdofinv,double ran);
/* ---- what did deleting those elements release (release.c) -------------

   Five sites in nonlingeo() and twelve locals owned by none of them, two of
   the sites being ninety-one-line arm blocks that were character-identical
   and sat 1,120 lines apart.  It reads and prints; it changes no equation
   and nothing downstream branches on it. */
typedef struct{
  ITG    probe;                /* 0 off, 1 measure, 2 per-element lines   */
  ITG    armed,pass;
  ITG    rebuild,iforbou;      /* what the topology event did             */
  ITG    nterm,nother,nisl,ncoh;   /* how the batch was made up           */
  double qa,qam,dt;            /* the reference force when it was armed   */
  double *frel;                /* f_int(u*) in NODE space, pre-remastruct */
  ITG    *ract;                /* which of those dofs were active then    */
}release;

void release_init(release *r);
ITG  release_configure(release *r);
void release_free(release *r);
void release_arm(release *r,const trialctx *mdl,const nlstate *n,
                 const erosion_batch *b,const topo_txn *txn,
                 const erosion_policy *pol,const double *damvisc,
                 ITG rebuild,ITG iforbou,ITG nisl,ITG ncoh);
void release_report(release *r,const trialctx *mdl,double steptime);

/* ---- where the Newton solve IS (nlstate.c) ----------------------------

   The iteration counter and the convergence quantities, in one object,
   because four different mechanisms ask for them and none of them owned
   them.  Every field is the ADDRESS of a local of nonlingeo(): the counter
   changes every iteration and the norms are rewritten in place, so a
   context holding values would be stale before its first use.  For the
   fixed-size arrays the array IS the address, so one rule covers both.

   Deliberately holds only what its readers ask for.  A context that tries
   to anticipate its callers becomes a second copy of nonlingeo()'s
   prologue, which is the thing being fixed.  Adding a field is one line
   here and one in NLSTATE_BIND. */
struct nlstate{
  ITG    *iit;                 /* Newton iteration inside this attempt   */
  double *ram,*ram1,*ram2;     /* residual norms: now, and two back      */
  double *cam;                 /* correction norms                       */
  double *uam;                 /* largest displacement increment         */
  double *qa;                  /* force quantities this iteration        */
  double *qam;                 /* the reference force                    */
  double *ctrl;                /* the *CONTROLS table                    */
};

/* Bind to the caller's frame.  Valid from the point of the call to the end
   of the function, no matter what is reallocated in between, because
   nothing here is a copy. */
#define NLSTATE_BIND(N) do{                                             \
  memset(&(N),0,sizeof(N));                                             \
  (N).iit=&iit;   (N).ram=ram;   (N).ram1=ram1; (N).ram2=ram2;          \
  (N).cam=cam;    (N).uam=uam;   (N).qa=qa;     (N).qam=qam;            \
  (N).ctrl=ctrl;                                                        \
}while(0)

ITG nlstate_check(const nlstate *n);

/* ---- which globalization mechanism did anything (globalize.c) ---------

   Six mechanisms stacked in a fixed order
   with no interface; it was recorded that three consecutive
   attempts producing bit-identical residual sequences, two rescue levels
   having run and changed nothing.  This turns that anecdote into a
   counter, so that "name the failure it addresses" can be answered with
   data before anything is deleted.

   Reason vocabulary from PETSc SNESLineSearchReason (petscsnes.h:918).
   It counts and reports; it changes no arithmetic.                     */

typedef enum {
  GLOB_LADDER      =0,
  GLOB_BACKTRACK   =1,
  GLOB_RESCUE1     =2,
  GLOB_RESCUE2     =3,
  GLOB_TRUSTREGION =4,
  GLOB_PATHFOLLOW  =5,
  GLOB_NMECH       =6
}glob_mech;

typedef struct{
  ITG fired[GLOB_NMECH];      /* the mechanism said it acted            */
  ITG effective[GLOB_NMECH];  /* ...and the correction changed          */
  ITG inert[GLOB_NMECH];      /* ...and it did not                      */
  ITG unknown[GLOB_NMECH];    /* ...on an attempt with no predecessor   */
  ITG attempts,bare_attempts,repeats;
  unsigned long long prev_hash,cur_hash;
  ITG have_prev,niter;
  unsigned pending;           /* mechanisms that acted this attempt     */
}glob_census;

const char *glob_mech_name(glob_mech m);
void glob_census_init(glob_census *g);
void glob_fired(glob_census *g,glob_mech m);
void glob_iterate(glob_census *g,const double *b,ITG n);
void glob_attempt_end(glob_census *g);
void glob_census_arm(const glob_census *g);
void glob_census_report(const glob_census *g);
ITG  glob_selftest(void);

/* ---- the transaction that commits an erosion (topology.c) -------------

   The marked set is one object with one
   lifetime, on the DMLabel pattern (PETSc include/petscdmlabel.h) rather
   than nine locals whose freeing was copy-pasted to four places.  It
   records a deletion the damage model has already decided; it decides
   nothing and it is not a rollback.                                    */

struct topo_txn{
  ITG *elem,*mat,*ip;      /* one entry per element eroded in this batch */
  double *value;
  ITG count;
  ITG step,increment;      /* when this batch happened                  */
  double step_time,total_time;
};

void topo_txn_init(topo_txn *t);
void topo_txn_discard(topo_txn *t);
ITG  topo_element_nip(const char *lakonel,ITG mi0);
void topo_txn_collect(topo_txn *t,ITG step,ITG increment,
                      double step_time,double total_time,
                      ITG ne0,const ITG *ipkondamageini,const ITG *ipkon,
                      const ITG *ielmat,const ITG *mi,const char *lakon,
                      const double *dam,
                      const ITG *ndmcon,const double *dmcon,
                      ITG ndmat_,ITG ntmat_,
                      ITG de13_transaction,const double *de13_trigger_value,
                      const ITG *de13_trigger_ip);
void topo_txn_write_history(const topo_txn *t,FILE *fdamage,ITG batch);
void topo_txn_commit(const topo_txn *t,FILE *fdamage,ITG batch,
                     ITG de13_transaction,ITG active_pass);
ITG  topo_selftest(void);

/* ---- who drives the load parameter (loadctl.c) ------------------------

   Four mechanisms that all answer the same question - what sets lambda
   this increment - and that all refuse to run beside each other:
   dissipation control, path control, the arc-length boundary, and the
   regularisation ladder.  Fifty-three locals of nonlingeo() with no owner,
   and the reason they matter beyond their own code is that every OTHER
   mechanism's arming block reads them: the dogleg and the continuation
   both refuse to arm when any of these is driving.  Those refusals used to
   name four loose integers; now they name one object.

   Fields keep their mechanism's prefix.  This is one object because the
   mutual exclusion is one rule, not because the four are one algorithm.

   `path_arm' was declared here and never read; it is gone.            */

typedef struct{
  /* the arc-length boundary: armed, and where it is */
  ITG    arc;
  double arc_lam,arc_theta0;
  /* dissipation control (CCX_DISSIPATION_*) */
  ITG    diss_ctrl,diss_engaged,diss_have,diss_init,diss_ok,diss_probe;
  ITG    diss_report,diss_step;
  char  *diss_env;
  double diss_den,diss_dg,diss_dgcur,diss_dgold,diss_dlam,diss_dtheta;
  double diss_engage_t,diss_ff,diss_fr,diss_g,diss_kpp,diss_lamcur;
  double diss_lamnow,diss_lamold,diss_lprev,diss_p,diss_pprev,diss_scale;
  double diss_slope,diss_target,diss_total;
  double *diss_fhat,*diss_uf;
  /* path control (CCX_DAMAGE_PATH): lambda as an absolute step fraction */
  ITG    path_att,path_desc,path_nstep,path_on,path_retry,path_used;
  double path_dev,path_devmax,path_drop,path_lam,path_lamcom,path_ref;
  /* the regularisation ladder */
  ITG    reg_level,reg_napply,reg_nlam,reg_on;
  double reg_lam[5],reg_lambda;
}loadctl;

void loadctl_init(loadctl *c);
/* Is any of the four driving the load parameter?  The one question the
   other mechanisms' arming blocks ask of this object. */
ITG  loadctl_driving(const loadctl *c);

/* ---- the operator check's driver state (opcheck.c) --------------------

   opcheck.c owns the comparison - the assembled tangent, column by column,
   against a central difference of the internal force - and the
   classification of what a discrepancy means.  What it did not own was the
   thirty-nine locals of nonlingeo() that arm it, hold the perturbed states
   it differentiates around, and tally the rank-1 census that says WHY the
   assembled tangent is not the differential of the residual.          */

typedef struct{
  /* the central-difference operator check (CCX_STRUCT_FD_*) */
  double  *fd_ad,*fd_au,*fd_f0,*fd_fm,*fd_fp,*fd_vsav,*fd_vtrue;
  ITG     fd_base,fd_col,fd_el,fd_inc,fd_it,fd_j,fd_ncol,fd_node,fd_pick,
          fd_s,fd_step,fd_t,fd_tel,fd_tnn,fd_uel;
  double  fd_dmax,fd_h,fd_udmax;
  /* the asymmetric tangent and its census (CCX_DAMAGE_TANGENT) */
  ITG     unsym_active,unsym_adv,unsym_advrep,unsym_census,unsym_degen,
          unsym_elems,unsym_floor,unsym_hole,unsym_holerep,unsym_live,
          unsym_report,unsym_skip,unsym_tanfull;   /* unsym_skiprep
          was declared in nonlingeo() and never read; it is gone */
}opcheckdrv;
void opcheckdrv_init(opcheckdrv *p);
void opcheckdrv_configure_fd(opcheckdrv *p);

/* ---- the probes' driver state (damdiag.c) -----------------------------

   damdiag.c owns ten pure observers.  This is what nonlingeo() holds
   between calls to them: which probe is armed, at which increment, the
   category maps they compare against, and the counters they fill.  Fifty
   locals, six clusters, one object - each field keeps its cluster's name
   because more than one cluster had an `inc'.                         */

typedef struct{
  /* the residual ray along the Newton direction */
  ITG     *ray_cat;
  double  ray_growth;
  ITG     ray_inc[4],ray_incok,ray_max,ray_ninc,ray_probe,ray_shots;
  double  *ray_p,*ray_r0,*ray_res;
  /* the wall probes: where the correction lives */
  ITG     wall_armed,wall_maskstep,wall_nb[8],wall_null;
  ITG     *wall_cat;
  double  *wall_def;
  double  wall_theta;
  /* the A-B-A full-state purity test */
  double  aba_a;
  ITG     aba_done,aba_hit,aba_inc[4],aba_mode,aba_ninc;
  /* the tension/compression event census */
  ITG     evt_fe[8],evt_fp[8],evt_nstep,evt_nsw[8],evt_on;
  ITG     *evt_sgn;
  /* the per-node dump */
  ITG     dump_alive,dump_hit,dump_idx,dump_inc,dump_n,dump_nb,dump_node,
          dump_np,dump_nu;
  double  dump_dv;
  /* the null-vector probe */
  double  null_amax,null_nb,null_nn,null_nx;
  ITG     null_cnt,null_inc,null_it,null_nit,null_seed;
  double  *null_x;
}probedrv;
void probedrv_init(probedrv *p);
void probedrv_configure_aba(probedrv *p);

/* ---- what happens when an increment will not converge (rescue.c) -----

   Sixty-seven locals of nonlingeo(), across six clusters that are one
   mechanism: the damage line search and its probe, transactional
   backtracking, the same-load re-equilibration after a deletion, the
   bounded recovery corridor, and the rescue levels that order them.

   They are grouped, not merged.  Each field keeps the name of the cluster
   it came from, because four of those clusters had a `mode' and a flat
   rename would have silently made them one field - a change that compiles
   and is wrong.  That is the failure mode of this transformation, and the
   reason the field list is printed and read before the struct is written.

   The LADDER that orders these levels is still in nonlingeo(): it cuts
   increments, rolls back topology and ends steps.  This is its state.  */

/* The damage line search: a full Newton correction is kept whenever it
   contracts the residual, and only a genuine increase in an active
   softening trial buys one safeguarded secant correction. */
#define DAMAGE_LINESEARCH_GROWTH 1.10
#define DAMAGE_LINESEARCH_MIN 0.10
#define DAMAGE_LINESEARCH_MAX 0.80
#define DAMAGE_LINESEARCH_FALLBACK 0.50
#define DAMAGE_LINESEARCH_MAX_TRIALS 3

typedef struct{
  /* transactional backtracking (CCX_DAMAGE_BT_*) */
  double  *bt_dam,*bt_visc,*bt_xs;
  double  bt_floor,bt_growth,bt_r0,bt_ring[8];
  ITG     bt_mode,bt_nring,bt_ntrial,bt_window;
  /* the bounded recovery corridor (CCX_DAMAGE_CORR_*) */
  ITG     corr_clean,corr_dtn,corr_exit,corr_grace,corr_maxesc,corr_maxinc;
  ITG     corr_maxwall,corr_mode,corr_nfact,corr_ninc,corr_nint,corr_nsince;
  ITG     corr_nwallstab,corr_nwalltot,corr_on,corr_stableneed,corr_trial;
  ITG     corr_try,corr_tryevery;
  double  corr_dtref,corr_dtsum,corr_lam,corr_lamstable,corr_minfrac;
  double  corr_theta0;
  /* the damage line search */
  ITG     linesearch_active,linesearch_applied,linesearch_contracted;
  ITG     linesearch_mode,linesearch_nsoft,linesearch_trial;
  double  linesearch_dampednorm,linesearch_fullnorm,linesearch_maxdd;
  double  linesearch_oldnorm;
  char    *linesearch_env;
  double  *linesearch_step;
  /* the line-search probe */
  ITG     ls_bestused,ls_legacy,ls_probe,ls_trials;
  double  ls_min;
  /* same-load re-equilibration after a deletion */
  char    *reeq_scale_env;
  ITG     reeq_scale_mode;
  double  reeq_uam_actual[2],reeq_uam_floor,reeq_uam_peak[2],reeq_uam_ref[2];
  /* the bounded recovery window: a rescue counts as recovered only
     after rec_window clean increments, and after rec_maxunrec
     un-recovered ones the ladder disarms itself */
  ITG    rec_disarmed,rec_healthy,rec_maxunrec,rec_unrec;
  ITG    rec_used_in_inc,rec_window;
  /* the rescue levels themselves */
  ITG     rescue_bt_on,rescue_maxlevel,rescue_mode,rescue_nfired;
  ITG     rescue_nok,rescue_used;
  double  rescue_dtheta_last,rescue_dthetaref_last;
}rescue;

void rescue_init(rescue *r);
/* the ladder's handshake with checkconvergence.c; defined in rescue.c */
extern ITG ccx_rescue_active,ccx_rescue_arm,ccx_rescue_req;
void rescue_configure_backtrack(rescue *r);
/* it WRITES to both: the corridor sets the regularisation ladder's
   length, and two levels reset the event census. */
void rescue_configure_levels(rescue *r,loadctl *c,probedrv *p);

/* ---- the Newton iteration budget (slownewton.c) -----------------------

   A decision that was two file-statics and thirteen locals: how many
   iterations may this increment have before the stock controller cuts it?
   The extension is bounded, and it applies only while erosion_softening()
   says damage is actually advancing.  It extends a budget; it relaxes no
   tolerance and closes no divergence path.                            */

#define DAMAGE_SLOW_NEWTON_MAX_ITERS 40
#define DAMAGE_SLOW_NEWTON_MAX_EXTRA 20
#define DAMAGE_SLOW_NEWTON_INVALID_EST 1000000000

typedef struct{
  ITG    active,allow,extended,maxiters;
  ITG    estres,estcorr,esttotal,nsoft;
  double camprev1,camprev2,cratio,rratio,maxdd;
}slownewton;

void slownewton_init(slownewton *s);
void slownewton_init(slownewton *s);
ITG slownewton_estimate(ITG iit,double value,
                                       double previous,double target);
/* Seventeen arguments became two.  Nine of them were the Newton state,
   which nlstate now owns; the other eight - camprev1, camprev2, maxiters
   and the five results - were ALREADY fields of this same slownewton
   object, taken apart at the call site and passed back in one at a time. */
ITG slownewton_allow(slownewton *s,const nlstate *n);

/* ---- what the damage state looks like, reported (damstats.c) ---------

   Three writers and the eight counters they fill.  None of them decides
   anything and the solver never reads back what they wrote.           */

/* DE1.1 fixed-point closure: the tolerance an accepted same-load pass has
   to meet, the relaxed one allowed at the pass limit, that limit, and the
   factor a failed closure cuts the physical increment by. */
#define DAMAGE_DE1_FP_TOL 1.e-4
#define DAMAGE_DE1_FP_RELAX_TOL 5.e-4
#define DAMAGE_DE1_MAX_PASSES 8
#define DAMAGE_DE1_CUTBACK_FACTOR 0.50

typedef struct{
  ITG    nactive,gt01,gt05,gt09,nfull,nchanged;
  double dmax,maxdelta;
}damstats;

void damstats_init(damstats *d);

/* Exact DE1 element statistics.
   For each element the maximum degradation over its active integration
   points is used.  In the present DE1 implementation only C3D4 is enabled,
   therefore this is exactly the single integration-point value. */
/* Fourteen arguments became three.  Five were trialctx's; the other eight
   were the fields of the damstats object the call was ABOUT, taken apart at
   the call site and handed back one at a time - the same shape the call to
   slownewton_allow() had.  damold stays a parameter because the three call
   sites pass different baselines. */
void damstats_element(damstats *d,const trialctx *mdl,const double *damold);
void damstats_append(const char *jobnamec,ITG istep,ITG iinc,
                                    double steptime,double totaltime,
                                    ITG passes,ITG nactive,ITG ngt01,
                                    ITG ngt05,ITG ngt09,ITG nfull,
                                    double dmax,double maxdelta);
/* Sixteen arguments became three: thirteen of them were the mesh, the
   material map and the damage state, all of which trialctx holds. */
void damstats_write_vtk(const char *jobnamec,const trialctx *mdl,
                        double steptime);

/* ---- the path-following DRIVER's state (declared here, driven from
   nonlingeo()) ---------------------------------------------------------

   pathfollow.c owns the method: the constraint, the bordered step, the
   predictor, the commit.  What it did not own was the sixty-three locals
   of nonlingeo() that drive it - whether it is armed, whether it has
   engaged, the frozen snapshot a rejected attempt is restored from, the
   nine workspace vectors, the linearity check and the mixed-mode crack
   control that rides on the same machinery.  Those are one object with
   one lifetime and they all already carried the prefix.

   The driver LOOP stays in nonlingeo(): it sets boundary conditions,
   decides cutbacks and ends increments.  This is its state, not its
   algorithm.                                                          */

typedef struct{
  /* armed, engaged, and why */
  ITG    on,engaged,applied,pending,reason,neqarm,lincheck,codmode;
  ITG    ncut,ncutfloor,nstep,icutbprev;
  char  *env;
  /* the constraint and the load factor it drives */
  double tauv,taucur,clip,eps,dtheta_eng;
  double lam,lamprev,lam0it,dlam,dlamit,dlamjump,dlampred;
  double g,dg,dgc,pdu,cu;
  double dphi,dphicur,gacc,phiacc;
  double fhcos,fhrat;
  double *cvec,*lhs,*rhs,*rhs0,*p,*q,*r0,*r1,*uf,*uref,*y;
  /* the snapshot a rejected attempt is restored from */
  double *sv,*sxs,*sxst,*sf,*sfn,*sstx,*sdam,*sxb;
  double sqa[4],scam[5];
  /* mixed-mode crack control, which rides on the same machinery */
  ITG    ccmode,ccnw,ccengage,ccarmed;
  double ccgrow;
  crackcontrol_census cs;
}pathdrv;

void pathdrv_init(pathdrv *p);
/* the driver's switches, read beside the method they drive; it refuses
   rather than degrades.  Declared after trialctx exists. */
double pathfollow_dot(const double *a,const double *b,ITG n);
double pathfollow_project(const double *fh,const double *a,const double *c,
                          const ITG *nactdof,ITG nk,ITG mt);

/* ---- the bounded local continuation (damcont.c) -----------------------

   Eighty-eight locals of nonlingeo() with no owner, across seventeen sites
   and 13,661 lines - the largest single cluster this fork added - are one
   object with one lifetime.  What was frozen, what the bordered system is
   solving, what the budgets are, what has been spent.

   The arithmetic (the bordered solve, the curvature denominator, the
   control-point choice) is here and has a self test against a closed form.
   The DRIVER - arming, the corrector loop, the commit lifecycle, the
   partial exit - is still in nonlingeo(); damcont.c says why.          */

typedef struct{
  /* Twelve more were declared here and never read - `have', `pred',
     `nstep', `nkink', `nretry', `lastminc', `taupr', `p1'..`p3' and two
     budgets.  They are gone.  Nothing in a function with a thousand
     locals could have shown that; an object can. */
  /* is it on, and how it was armed */
  ITG    on,arm,mode,alloc,used,refused,partial;
  /* the frozen control point and its kinematics */
  ITG    elem,ip,head,nring,nsupp,nbl,p1r,p2r,p3r;
  ITG    node[18],dir[18],supp[18];
  double m[3],g[3],dc[3];
  double duinf[6],dt[6],dlam[6];
  double clam,cprev,duref,lam,lamc,lamref,lamsnap;
  double qa[4],cam[5],uam[2];
  /* the bordered system and the step */
  double ds,ds0,dsmin,dsmax,den,rho,tau,kappa,tolc;
  double eps,kaptol,rhomin,clim,ulim;
  ITG    epsok;
  double *beps,*dam,*jac,*qh,*r0,*ring,*visc,*w,*xs,*y,*z;
  ITG    *bl,*fl,*sgn;
  /* budgets, and what has been spent against them */
  ITG    maxcorr,maxeval,maxfact,maxstep;
  ITG    it,step,ncorr,neval,nfact,nprog,ncommit;
  ITG    newstep;
}damcont;

void damcont_init(damcont *c);
void damcont_init(damcont *c);
void damcont_kin(const trialctx *mdl,ITG indexe,const double *v,ITG mint,
                 double *dl,double *rmat,double *shape);
/* Thirteen arguments became five.  v and stx stay parameters: which state
   the caller wants snapped is its decision, and the one call site passes
   vold and sti rather than the context's v and stx. */
void damcont_snap(const trialctx *mdl,const double *v,const double *stx,
                  double *ring,ITG *fl);
ITG damcont_bordered(double clam,double cuz,double cuy,double c,
                              double *den,double *dlam);
double damcont_rhoden(double clam,double cunorm,double ysupp,
                               double den);
ITG damcont_select(const double *ring,const ITG *fl,
                            const double *dtr,const ITG *ipkon,
                            const char *lakon,const ITG *ielprop,
                            const double *prop,ITG ne0,ITG mi0,ITG head,
                            ITG *belem,ITG *bip,double *m,double *kappa,
                            double *ds0,double *tau,ITG *ncand,ITG *reason,
                            double kaptol,const ITG *bl,ITG nbl);
ITG damcont_selftest(void);

/* ---- the trust region and its state (dogleg.c) ------------------------

   Forty-nine locals of nonlingeo() with no owner, across twenty sites and
   15,498 lines, are one object: they are created together, they mean
   nothing apart, and every one of them was already called damage_dl_*.
   Policy, workspace and census, in that order.

   dogleg_pick() is the step CHOICE alone - a radius and five scalars in,
   a branch and two coefficients out - which is what makes it testable
   against a closed form.  The LOOP that accepts, rejects and resizes is
   still in nonlingeo(); see the block comment in dogleg.c for why, and
   for what changed about whether it has to be.                        */

typedef struct{
  /* policy - read from the switches at arming */
  ITG    on,mode,banner,lincheck,selfrec,incarm;
  ITG    maxarm,maxeval,maxfact,maxtrial;
  double d0fac;
  /* workspace - allocated when the region arms, freed with it */
  double *d,*pn,*pm,*r0,*res,*w,*wm;
  double *xs,*dam,*visc;          /* the snapshot a trial is rolled back to */
  /* the current step: the five scalars dogleg_pick() works in, and the
     radius and residual they belong to */
  double nb2,nd2,nw2,npn2,npm2,dtpn,tc,delta,phi0,dmax,ident,asym;
  ITG    have;
  /* census - what the mechanism actually did, over the run */
  ITG    narm,nacc,nrej,nnewt,ncau,ndog,neval,nfact,nfail,used;
  ITG    recdone,lasthelp,lc_due,lc_it,lc_nit;
}dogleg;

void dogleg_init(dogleg *d);
/* the step: 1 NEWTON, 2 CAUCHY, 3 DOGLEG, 0 REFUSE (degenerate or NaN -
   which the caller must treat as no step, never as a step of length 0) */
ITG  dogleg_pick(double dl,double nd2,double nw2,double npn2,
                 double dtpn,double *pa,double *pb,double *nrm);
/* the model's predicted reduction for that step: the denominator of rho */
double dogleg_pred(double pa,double pb,double nb2,double nd2,double nw2);
ITG  dogleg_selftest(void);
/* read this mechanism's own switches, and refuse rather than degrade.
   Called where the block used to be: these arming blocks refuse on
   each other's state, so their order is part of the behaviour. */
void dogleg_configure(dogleg *d,rescue *r,const loadctl *c);
/* level 4 arms on top of the dogleg, so its declaration lives here,
   after both types exist. */
void damcont_configure(damcont *k,dogleg *d,rescue *r,
                       const loadctl *c);

/* ---- the probes (damdiag.c) -------------------------------------------

   Ten pure observers that were file-statics of nonlingeo.c.  They print
   and they return numbers; they change no solver state, which is the
   whole boundary and the reason moving them could not change an answer.
   damage_ray_catof() stays private to that file: what a category IS
   belongs with the three functions that ask.                          */

/* ---- discrete-branch census -------------------------------------------
   J-14 measured a sharp loss of local linearity between alpha=0.0625 and
   alpha=0.125 and called it a "discrete switch".  That was a hypothesis, not
   a measurement: nothing had been counted.  These two helpers count it.

   One bitmask per integration point.  The element type decides which bits are
   meaningful, and getting that wrong is not hypothetical - xstate slot 1 is
   the equivalent plastic strain for a bulk C3D4 (calcdamage.f:649) and dmax
   for a UC6 facet (resultsmech_uc6.f:43).  Reading one as the other yields a
   plausible, meaningless census.

   Layout: dam(mi(1),*)          -> dam[mi0*i+j]
           xstate(nstate_,mi(1),*) -> xstate[nstate_*mi0*i + nstate_*j + k] */

#define DAMCAT_PLAST   1   /* bulk: accumulated plastic strain this increment */
#define DAMCAT_DINIT   2   /* bulk: damage initiated (dam >= 1)               */
#define DAMCAT_DGROW   4   /* bulk: damage grew from the committed baseline   */
#define DAMCAT_USOFT   8   /* UC6 : dmax > 0, law has left the elastic branch */
#define DAMCAT_UVISC  16   /* UC6 : viscous damage active                     */
#define DAMCAT_UFAIL  32   /* UC6 : fully failed flag set                     */
#define DAMCAT_UADV   64   /* UC6 : dmax ADVANCED this trial, i.e. deff>dmax0  */
#define DAMCAT_UCOMP 128   /* UC6 : traction(1)<0, the compression branch     */

void damage_aba_cmp(const char *name,const double *a,
                    const double *b,ITG n,ITG *nbad);
/* Four functions in this file shared an identical eleven-parameter prefix -
   the state arrays, the mesh and the sizes.  All of it is trialctx's except
   dambase and visc. */
void damage_ray_census(ITG *cat,const trialctx *mdl,
                       const double *dambase,const double *visc);
void damage_ray_tally(const trialctx *mdl,const double *dambase,
                      const double *visc,ITG *nplast,ITG *nucomp);
void damage_evt_sign(const trialctx *mdl,ITG *sgn);
ITG damage_evt_flips(const trialctx *mdl,const ITG *sgn,
                     ITG *firste,ITG *firstip);
ITG damage_ray_census_diff(const ITG *cat,const trialctx *mdl,
                           const double *dambase,const double *visc,
                           ITG *firste,ITG *firstip,
                           ITG *firsta,ITG *firstb);
/* Seventeen arguments became four: thirteen were the mesh and the state,
   which trialctx holds, and mt is mi[1]+1. */
void damage_wall_where(const char *tag,const double *x,const trialctx *mdl,
                       ITG ntop);
/* Fifteen arguments became three, for the same reason. */
void damage_wall_split(const char *tag,const double *x,const trialctx *mdl);
ITG damage_wall_setdiff(const ITG *cat,const trialctx *mdl,
                        const double *dambase,const double *visc,ITG *nb);

/* ---- evaluating the residual at a trial state (trial.c) ---------------

   The one operation every globalisation mechanism here is built from, and
   the one that was written out by hand seventeen times inside nonlingeo():
   put a step in b, evaluate the model, read the residual.

   trialctx is NOT an abstraction of the model.  It is the argument list of
   results() and calcresidual() written down once - one field per argument,
   each holding the ADDRESS of the caller's local so that the binding
   survives every NNEW, SFREE and remastruct.  The struct, TRIAL_BIND and
   the calls in trial.c were GENERATED from the call text and these
   prototypes; the rule is uniform (`&x' for an argument written `x' or
   `&x') which is what makes trial_check() able to catch a field the bind
   forgot.  See the block comment at the top of trial.c.              */

struct trialctx{
  double **    co;
  ITG **       nk;
  ITG **       kon;
  ITG **       ipkon;
  char **      lakon;
  ITG **       ne;
  double **    v;
  double **    stn;
  ITG **       inum;
  double **    stx;
  double **    elcon;
  ITG **       nelcon;
  double **    rhcon;
  ITG **       nrhcon;
  double **    alcon;
  ITG **       nalcon;
  double **    alzero;
  ITG **       ielmat;
  ITG **       ielorien;
  ITG **       norien;
  double **    orab;
  ITG **       ntmat_;
  double **    t0;
  double **    t1act;
  ITG **       ithermal;
  double **    prestr;
  ITG **       iprestr;
  char **      filab;
  double **    eme;
  double **    emn;
  double **    een;
  ITG **       iperturb;
  double **    f;
  double **    fn;
  ITG **       nactdof;
  ITG *        iout;
  double *     qa;   /* an array local: it never moves */
  double **    vold;
  double **    b;
  ITG **       nodeboun;
  ITG **       ndirboun;
  double **    xbounact;
  ITG **       nboun;
  ITG **       ipompc;
  ITG **       nodempc;
  double **    coefmpc;
  char **      labmpc;
  ITG **       nmpc;
  ITG **       nmethod;
  double *     cam;   /* an array local: it never moves */
  ITG *        neq1;
  double **    veold;
  double **    accold;
  double *     bet;
  double *     gam;
  double *     dtime;
  double *     time;
  double **    ttime;
  double **    plicon;
  ITG **       nplicon;
  double **    plkcon;
  ITG **       nplkcon;
  double **    xstateini;
  double **    xstiff;
  double **    xstate;
  ITG **       npmat_;
  double **    epn;
  char **      matname;
  ITG **       mi;
  ITG *        ielas;
  ITG *        icmd;
  ITG **       ncmat_;
  ITG **       nstate_;
  double **    stiini;
  double **    vini;
  ITG **       ikboun;
  ITG **       ilboun;
  double **    ener;
  double **    enern;
  double **    emeini;
  double **    xstaten;
  double **    eei;
  double **    enerini;
  double **    cocon;
  ITG **       ncocon;
  char **      set;
  ITG **       nset;
  ITG **       istartset;
  ITG **       iendset;
  ITG **       ialset;
  ITG **       nprint;
  char **      prlab;
  char **      prset;
  double **    qfx;
  double **    qfn;
  double **    trab;
  ITG **       inotr;
  ITG **       ntrans;
  double **    fmpc;
  ITG **       nelemload;
  ITG **       nload;
  ITG **       ikmpc;
  ITG **       ilmpc;
  ITG **       istep;
  ITG *        iinc;
  double **    springarea;
  double *     reltime;
  ITG *        ne0;
  double **    thicke;
  double **    shcon;
  ITG **       nshcon;
  char **      sideload;
  double **    xloadact;
  double **    xloadold;
  ITG *        icfd;
  ITG **       inomat;
  double **    pslavsurf;
  double **    pmastsurf;
  ITG **       mortar;
  ITG **       islavact;
  double **    cdn;
  ITG **       islavnode;
  ITG **       nslavnode;
  ITG **       ntie;
  double **    clearini;
  ITG **       islavsurf;
  ITG **       ielprop;
  double **    prop;
  double *     energyini;   /* an array local: it never moves */
  double **    energy;
  ITG *        kscale;
  ITG **       iponoeln;
  ITG **       inoeln;
  ITG **       nener;
  char **      orname;
  ITG **       network;
  ITG **       ipobody;
  double **    xbodyact;
  ITG **       ibody;
  char **      typeboun;
  ITG **       itiefac;
  char **      tieset;
  double **    smscale;
  ITG *        mscalmethod;
  ITG **       nbody;
  double **    t0g;
  double **    t1g;
  ITG **       islavquadel;
  double **    aut;
  ITG **       irowt;
  ITG **       jqt;
  ITG *        mortartrafoflag;
  ITG *        intscheme;
  double **    physcon;
  double **    dam;
  double **    damn;
  ITG **       iponoel;
  ITG **       neq;
  double **    res;
  double **    fext;
  ITG **       iexpl;
  double **    aux2;
  double **    adb;
  double **    aub;
  ITG **       jq;
  ITG **       irow;
  ITG **       nzl;
  double **    alpha;
  double **    fextini;
  double **    fini;
  ITG **       nzs;
  ITG *        nasym;
  ITG *        idamping;
  double **    adc;
  double **    auc;
  double **    cvini;
  double **    cv;
  double *     alpham;
  ITG *        num_cpus;
  /* not an argument of either call, but the atom brackets results()
     with it: the 1d/2d expansion needs inum and nothing else does. */
  ITG *        ne1d2d;
};

/* Bind the context to the caller's frame.  One line per field, next to the
   struct so the two cannot drift apart; trial_check() catches it if they
   do.  Call it once, after every local it names exists. */
/* Zero FIRST.  Without it the check below this macro cannot do the job its
   own comment claims: `trialctx nlgt;' is an uninitialised local, so a field
   TRIAL_BIND forgot holds stack garbage, not NULL, and trial_check() walks
   past it reporting 0 unbound.  OBSERVED - commenting a field out of the
   bind left the check green.  One memset makes the claim true. */
#define TRIAL_BIND(T) do{ \
  memset(&(T),0,sizeof(T)); \
  (T).co=&co; \
  (T).nk=&nk; \
  (T).kon=&kon; \
  (T).ipkon=&ipkon; \
  (T).lakon=&lakon; \
  (T).ne=&ne; \
  (T).v=&v; \
  (T).stn=&stn; \
  (T).inum=&inum; \
  (T).stx=&stx; \
  (T).elcon=&elcon; \
  (T).nelcon=&nelcon; \
  (T).rhcon=&rhcon; \
  (T).nrhcon=&nrhcon; \
  (T).alcon=&alcon; \
  (T).nalcon=&nalcon; \
  (T).alzero=&alzero; \
  (T).ielmat=&ielmat; \
  (T).ielorien=&ielorien; \
  (T).norien=&norien; \
  (T).orab=&orab; \
  (T).ntmat_=&ntmat_; \
  (T).t0=&t0; \
  (T).t1act=&t1act; \
  (T).ithermal=&ithermal; \
  (T).prestr=&prestr; \
  (T).iprestr=&iprestr; \
  (T).filab=&filab; \
  (T).eme=&eme; \
  (T).emn=&emn; \
  (T).een=&een; \
  (T).iperturb=&iperturb; \
  (T).f=&f; \
  (T).fn=&fn; \
  (T).nactdof=&nactdof; \
  (T).iout=&iout; \
  (T).qa=qa; \
  (T).vold=&vold; \
  (T).b=&b; \
  (T).nodeboun=&nodeboun; \
  (T).ndirboun=&ndirboun; \
  (T).xbounact=&xbounact; \
  (T).nboun=&nboun; \
  (T).ipompc=&ipompc; \
  (T).nodempc=&nodempc; \
  (T).coefmpc=&coefmpc; \
  (T).labmpc=&labmpc; \
  (T).nmpc=&nmpc; \
  (T).nmethod=&nmethod; \
  (T).cam=cam; \
  (T).neq1=&neq[1]; \
  (T).veold=&veold; \
  (T).accold=&accold; \
  (T).bet=&bet; \
  (T).gam=&gam; \
  (T).dtime=&dtime; \
  (T).time=&time; \
  (T).ttime=&ttime; \
  (T).plicon=&plicon; \
  (T).nplicon=&nplicon; \
  (T).plkcon=&plkcon; \
  (T).nplkcon=&nplkcon; \
  (T).xstateini=&xstateini; \
  (T).xstiff=&xstiff; \
  (T).xstate=&xstate; \
  (T).npmat_=&npmat_; \
  (T).epn=&epn; \
  (T).matname=&matname; \
  (T).mi=&mi; \
  (T).ielas=&ielas; \
  (T).icmd=&icmd; \
  (T).ncmat_=&ncmat_; \
  (T).nstate_=&nstate_; \
  (T).stiini=&stiini; \
  (T).vini=&vini; \
  (T).ikboun=&ikboun; \
  (T).ilboun=&ilboun; \
  (T).ener=&ener; \
  (T).enern=&enern; \
  (T).emeini=&emeini; \
  (T).xstaten=&xstaten; \
  (T).eei=&eei; \
  (T).enerini=&enerini; \
  (T).cocon=&cocon; \
  (T).ncocon=&ncocon; \
  (T).set=&set; \
  (T).nset=&nset; \
  (T).istartset=&istartset; \
  (T).iendset=&iendset; \
  (T).ialset=&ialset; \
  (T).nprint=&nprint; \
  (T).prlab=&prlab; \
  (T).prset=&prset; \
  (T).qfx=&qfx; \
  (T).qfn=&qfn; \
  (T).trab=&trab; \
  (T).inotr=&inotr; \
  (T).ntrans=&ntrans; \
  (T).fmpc=&fmpc; \
  (T).nelemload=&nelemload; \
  (T).nload=&nload; \
  (T).ikmpc=&ikmpc; \
  (T).ilmpc=&ilmpc; \
  (T).istep=&istep; \
  (T).iinc=&iinc; \
  (T).springarea=&springarea; \
  (T).reltime=&reltime; \
  (T).ne0=&ne0; \
  (T).thicke=&thicke; \
  (T).shcon=&shcon; \
  (T).nshcon=&nshcon; \
  (T).sideload=&sideload; \
  (T).xloadact=&xloadact; \
  (T).xloadold=&xloadold; \
  (T).icfd=&icfd; \
  (T).inomat=&inomat; \
  (T).pslavsurf=&pslavsurf; \
  (T).pmastsurf=&pmastsurf; \
  (T).mortar=&mortar; \
  (T).islavact=&islavact; \
  (T).cdn=&cdn; \
  (T).islavnode=&islavnode; \
  (T).nslavnode=&nslavnode; \
  (T).ntie=&ntie; \
  (T).clearini=&clearini; \
  (T).islavsurf=&islavsurf; \
  (T).ielprop=&ielprop; \
  (T).prop=&prop; \
  (T).energyini=energyini; \
  (T).energy=&energy; \
  (T).kscale=&kscale; \
  (T).iponoeln=&iponoeln; \
  (T).inoeln=&inoeln; \
  (T).nener=&nener; \
  (T).orname=&orname; \
  (T).network=&network; \
  (T).ipobody=&ipobody; \
  (T).xbodyact=&xbodyact; \
  (T).ibody=&ibody; \
  (T).typeboun=&typeboun; \
  (T).itiefac=&itiefac; \
  (T).tieset=&tieset; \
  (T).smscale=&smscale; \
  (T).mscalmethod=&mscalmethod; \
  (T).nbody=&nbody; \
  (T).t0g=&t0g; \
  (T).t1g=&t1g; \
  (T).islavquadel=&islavquadel; \
  (T).aut=&aut; \
  (T).irowt=&irowt; \
  (T).jqt=&jqt; \
  (T).mortartrafoflag=&mortartrafoflag; \
  (T).intscheme=&intscheme; \
  (T).physcon=&physcon; \
  (T).dam=&dam; \
  (T).damn=&damn; \
  (T).iponoel=&iponoel; \
  (T).neq=&neq; \
  (T).res=&res; \
  (T).fext=&fext; \
  (T).iexpl=&iexpl; \
  (T).aux2=&aux2; \
  (T).adb=&adb; \
  (T).aub=&aub; \
  (T).jq=&jq; \
  (T).irow=&irow; \
  (T).nzl=&nzl; \
  (T).alpha=&alpha; \
  (T).fextini=&fextini; \
  (T).fini=&fini; \
  (T).nzs=&nzs; \
  (T).nasym=&nasym; \
  (T).idamping=&idamping; \
  (T).adc=&adc; \
  (T).auc=&auc; \
  (T).cvini=&cvini; \
  (T).cv=&cv; \
  (T).alpham=&alpham; \
  (T).num_cpus=&num_cpus; \
  (T).ne1d2d=&ne1d2d; \
  }while(0)

/* WHICH OF THESE MOVES THE SCRATCH ARRAYS.  A caller that binds a local
   alias to *(mdl->v), *(mdl->stx), *(mdl->fn) or *(mdl->inum) needs to know, and
   the answer is: trial_results() and trial_reduce() do NOT - they read and
   write through the arrays that are there.  trial_evaluate() and
   trial_residual() free and reallocate all four, so an alias taken before
   one of them is dangling after it.  opcheck_probe() aliases v and fn and
   is correct because it only ever calls trial_results(); rescue_backtrack()
   reads stx through the context at each use because it does not have that
   luxury. */
void trial_results(const trialctx *mdl);          /* evaluate the model   */
void trial_evaluate(const trialctx *mdl);         /* ...with its scratch: REALLOCATES v,stx,fn,inum */
void trial_reduce(const trialctx *mdl,double *dst);/* ...reduce to a residual */
void trial_residual(const trialctx *mdl,double *dst);/* scratch + both halves: REALLOCATES */
ITG  trial_check(const trialctx *mdl);
/* the path follower's predictor: capture f_hat, solve for the reference
   direction, place lambda for this attempt.  It calls the linear solver,
   hence the factorisation arguments.  The caller keeps the guard. */
void pathdrv_predictor(pathdrv *p,glob_census *g,const trialctx *mdl,
                       const double *xboun,const double *xbounold,
                       double *ad,double *au,ITG *icol,const ITG *isolver,
                       double sigma,ITG inputformat,ITG nrhs,
                       ITG symmetryflag,ITG iit);
/* one Newton iteration of the bordered system, once level 4 owns the
   boundary.  The caller keeps the guard. */
void damcont_corrector(damcont *k,const trialctx *mdl,
                       const double *xboun,const double *xbounold,
                       double *uam,
                       const double *damjac,const double *damvisc,
                       ITG inputformat,ITG nrhs,ITG symmetryflag);
/* the one extra attempt a deferred stop buys, and which level gets it.
   Thirteen parameters: six objects and seven pieces of nonlingeo's own
   control flow.  See the block comment in rescue.c for why they are not
   bundled. */
void rescue_attempt(rescue *r,dogleg *d,damcont *k,
                    loadctl *c,probedrv *p,glob_census *g,
                    const trialctx *mdl,
                    double *dtheta,double *dthetaref,
                    double theta,const double *tper,ITG iit,ITG idamagereeq);
/* the operator check's measurement: the assembled tangent against a
   central difference, split by population.  A diagnostic that ends the
   run; the caller keeps the guard. */
void opcheck_probe(opcheckdrv *o,const trialctx *mdl,const ITG *ndmat_,
                   const ITG *damcat,ITG iit);
void pathdrv_configure(pathdrv *p,const loadctl *c,const trialctx *mdl,
                       const ITG *isolver,ITG ncont);

/* The trust-region LOOP.  239 lines that could not be moved before
   trialctx and the dogleg object existed, because there was no signature
   to give them; now there is one, of seven arguments.  The caller keeps
   the guard - whether the region may fire at all is a decision about the
   increment - and this does what the region then does. */
void dogleg_rescue(dogleg *d,const trialctx *mdl,glob_census *g,
                   double *damvisc,ITG iit,ITG icutb,double *uam);

/* Transactional backtracking, the second loop that could not be moved.
   Walks alpha down 1, 1/2 ... 1/64 with the committed baseline restored
   before every probe, accepts on Armijo against a non-monotone reference,
   restores the full step if nothing is acceptable.  It CHANGES THE ANSWER
   and is off by default; the caller keeps the guard. */
void rescue_backtrack(rescue *r,const trialctx *mdl,glob_census *g,
                      probedrv *p,double *damvisc,ITG iit);

/* ---- which elements leave the assembly (erosion.c) --------------------

   Three rules - the terminal damage threshold, DEADALL and DEADSOLE -
   applied in that order against one batch budget.  The block that applied
   them was written out twice in nonlingeo.c, at the two points where a
   converged state is scanned, identical but for indentation; erosion_mark()
   is the call both sites make now.  The three rules are file-static inside
   erosion.c, so "what else can delete an element" has the answer `nothing,
   by construction'.

   It MARKS (ipkon -> -ipkon-2) and does no more: the equation structure,
   the .damage record and the rollback belong to topology.c and to
   nonlingeo, and erosion.c contains no code that touches them.          */

struct erosion_policy{
  double delete_d;      /* CCX_DAMAGE_DELETE_D: the terminal threshold    */
  ITG    delete_visc;   /* CCX_DAMAGE_DELETE_VISC: read Dvis, not D       */
  const char *filter;   /* CCX_DAMAGE_DELETE_MAT; NULL means every one    */
  double deadall_g;     /* CCX_DAMAGE_DEADALL, 0 = the rule is off        */
  double deadsole_g;    /* CCX_DAMAGE_DEADSOLE, 0 = the rule is off       */
  ITG    batchmax;      /* how many may leave in one transaction          */
};

/* What one call took, and what the run has taken so far.  Seven locals of
   nonlingeo() with no owner, read from five other places; one object. */
struct erosion_batch{
  ITG marked;           /* elements marked by this call                   */
  ITG terminal;         /* ...of them by the damage threshold alone       */
  ITG deadall,deadsole; /* ...by each load-path rule                      */
  ITG deadall_nodes;    /* nodes DEADALL judged unsupported               */
  double batch_dmax;    /* worst D in the batch                           */
  double batch_vmin;    /* smallest Dvis in it - how much stress it held  */
  ITG total_deadall,total_deadsole;   /* running, over the whole run      */
};

void erosion_batch_init(erosion_batch *b);
ITG  erosion_mark(const erosion_policy *p,erosion_batch *b,
                  double *dam,const double *damvisc,
                  ITG *ipkon,const char *lakon,const ITG *kon,
                  const ITG *ielmat,const char *matname,
                  const ITG *ndmcon,const double *dmcon,
                  ITG ndmat,ITG ntmat,ITG nk,ITG ne,ITG ne0,ITG mi0,ITG mi2,
                  double *trigger_value,ITG *trigger_ip,
                  ITG iinc,double steptime);
/* Is anything actually softening in the present Newton trial?  The line
   search and the iteration budget both ask; having a DE1/DM2.0 material in
   the model is not an answer. */
/* The material table - ndmcon, dmcon, ndmat, ntmat - stays loose.  It is
   not trialctx's: trialctx transcribes results()'s argument list, and the
   damage constants are not among them. */
ITG  erosion_softening(const trialctx *mdl,const double *dambase,
                       const ITG *ndmcon,const double *dmcon,
                       ITG ndmat,ITG ntmat,
                       ITG *nsoft,double *maxdd);
ITG  erosion_selftest(void);

/* Whether a material's damage is progressive is a Material question and
   that object does not exist yet; erosion.c owns the answer meanwhile and
   topology.c calls it directly rather than inventing a dependency. */
ITG dammat_selftest(void);

/* The self test table, and the one place that runs it.  selftest_gate() is
   the interlock nonlingeo() calls: silent on a sound build, and it stops
   the run on a bad one.  ccx_selftest calls selftest_run_all(1) for the
   whole table with output. */
typedef ITG (*selftest_fn)(void);
typedef struct { const char *name; selftest_fn fn; } selftest_entry;
extern const selftest_entry SELFTESTS[];
extern const ITG NSELFTEST;
ITG  selftest_run_all(ITG verbose);
void selftest_gate(void);
ITG damage_progressive_material(ITG imat,const ITG *ndmcon,
                                const double *dmcon,ITG ndmat,ITG ntmat);

/* ---- what counts as converged (converge.c) ---------------------------

   The mechanical criterion is eight clauses joined by && and ||, written
   out three times in checkconvergence.c (mechanical, thermal, thermo-
   mechanical) with the names thrown away.  A NOX StatusTest Combo tree is
   exactly this shape with the names kept, so each leaf here gets a name, a
   value, the threshold it was compared against, and a status - which makes
   "which clause is holding this increment back" answerable.  The boolean
   is unchanged and bit-identity is the standard.                        */

typedef enum {
  CVG_UNEVALUATED=0,   /* the variant in force did not ask this leaf     */
  CVG_PASS       =1,
  CVG_FAIL       =2
}cvg_status;

typedef struct{
  const char *name;
  double value;        /* the quantity tested                            */
  double thresh;       /* what it was tested against                     */
  cvg_status status;
}cvg_clause;

#define CVG_MAXCLAUSE 20

typedef struct{
  ITG converged;              /* the boolean checkconvergence used to set */
  ITG nclause;
  cvg_clause clause[CVG_MAXCLAUSE];
  const char *blocker;        /* first failing leaf of the top-level AND  */
  double blocker_value,blocker_thresh;
}cvg_verdict;

/* the tolerances the verdict applies, named once instead of ctrl[] indices
   spelled out at each use site */
typedef struct{
  double ran,can,rap,ea,cae,ral,cetol;
  ITG ip;
}cvg_tol;

void cvg_tol_from_ctrl(cvg_tol *t,const double *ctrl);
ITG  converge_verdict(cvg_verdict *v,const cvg_tol *t,ITG ithermal,ITG iit,
                      ITG nmethod,ITG iflagact,ITG ntg,double deltmx,
                      const double *ram,const double *ram1,double *ram2,
                      const double *cam,const double *uam,
                      const double *qa,const double *qam,
                      double *c1,double *c2);
void converge_verdict_print(const cvg_verdict *v);

/* ---- why the run stopped (converge.c) --------------------------------

   PETSc's SNESConvergedReason exists because "it did not converge" is not
   an answer.  This tree ends with rc=201 and one message - "increment size
   smaller than minimum" - printed from four different places that mean
   three different things: the increment CONVERGED but wanted a next step
   below tmin; the residual DIVERGED and the cutback hit tmin; convergence
   was too SLOW and the cutback hit tmin.  Same exit code, same sentence.

   Sign convention follows PETSc: positive converged, negative diverged,
   zero still iterating.                                                */

typedef enum {
  CVG_DIVERGED_MINSTEP_EXTERNAL      =-4,
  CVG_DIVERGED_MINSTEP_TOO_SLOW      =-3,
  CVG_DIVERGED_MINSTEP_ON_DIVERGENCE =-2,
  CVG_DIVERGED_MINSTEP_AFTER_CONV    =-1,
  CVG_ITERATING                      = 0,
  CVG_CONVERGED_CRITERIA             = 1
}cvg_reason;

const char *converge_reason_name(cvg_reason r);
const char *converge_reason_explain(cvg_reason r);
void        converge_stop_report(cvg_reason r,const cvg_verdict *v);
ITG  converge_selftest(void);

/* ---- the backtracking ladder of the damage line search (lsladder.c) --
   Extracted from the Newton loop because it was wrong and the way it was
   wrong is worth a regression test.  See the block comment there. */

typedef struct{
  double alpha;     /* the trial step length now under test              */
  double best;      /* the best rung measured so far                     */
  double bestres;   /* and the residual there                            */
  double oldres;    /* the residual the search has to beat               */
  double floor;     /* the smallest rung the ladder may reach            */
  double ratio;     /* how fast the ladder descends                      */
  ITG ntrial;       /* rungs allowed                                     */
  ITG itrial;       /* rungs used                                        */
  ITG contracted;   /* a rung beat oldres                                */
  ITG legacy;       /* reproduce the old shape and the old fallback      */
}lsladder;

void lsladder_start(lsladder *l,double alpha0,double oldres,double flr,
                    double ratio,ITG ntrial,ITG legacy);
ITG lsladder_step(lsladder *l,double res);
double lsladder_final(const lsladder *l);
ITG lsladder_selftest(void);

/* ---- topology diagnostic (topodiag.c) -------------------------------
   Measurement only: nothing here decides anything or is on a solution
   path.  It exists to tell a true orphan degree of freedom, a physically
   detached component with rigid-body modes, and a sound topology with a
   failing corrector apart from one another. */

typedef struct{
  ITG nelem;        /* live elements the graph was built from            */
  ITG nfacet;       /* of those, user (cohesive) elements                */
  ITG ncomp;        /* connected components of the live element graph    */
  ITG maincomp;     /* id of the largest component                       */
  ITG mainnode;     /* nodes in it                                       */
  ITG nfloat;       /* components with NO prescribed dof and NO MPC      */
  ITG nfloatnode;   /* nodes in those                                    */
  ITG nfloatdof;    /* active dofs in those                              */
  ITG floatmax;     /* nodes in the largest floating component           */
  ITG floatfirst;   /* id of the first floating component                */
  ITG norphandof;   /* active dofs on nodes with no live element         */
  ITG norphannode;  /* nodes carrying them                               */
  ITG orphannode;   /* the first such node, 1-based, or -1               */
  ITG nzerodiag;    /* exactly zero diagonal entries                     */
  ITG ntinydiag;    /* below 1e-12 of the largest                        */
  ITG nisolated;    /* equations with no off-diagonal coupling at all    */
  ITG resmaxeq;     /* equation carrying the largest residual            */
  ITG resmaxnode;   /* and its node, 1-based                             */
  ITG resmaxdir;    /* and direction                                     */
  ITG resmaxcomp;   /* and component                                     */
  double admax,admin;
  double resnorm,resmax,resfloat;
}topodiag_report;

ITG topodiag_find(ITG *p,ITG a);
void topodiag_union(ITG *p,ITG a,ITG b);
ITG topodiag_nope(const char *lak);
ITG topodiag_selftest(void);
void topodiag_report_zero(topodiag_report *r);
/* Twenty-two arguments became five.  Seventeen of them were the mesh, the
   constraints and the equation structure, which trialctx has held all along;
   mt is mi[1]+1, which the caller was computing by hand.  Only the assembled
   diagonal and off-diagonal are still loose, and they are the linear system,
   which has no object yet. */
void topodiag_run(topodiag_report *r,ITG *comp,const trialctx *mdl,
                  const double *ad,const double *au);
double topodiag_project(const double *v,const double *w,ITG neq);
void topodiag_support(ITG node,const trialctx *mdl,ITG *nbulk,ITG *nfac);
ITG topodiag_deflate(double *w,const double *q,ITG k,ITG neq);
double topodiag_project_span(const double *v,const double *q,ITG k,
                             ITG neq);
void topodiag_print(const topodiag_report *r,const char *tag,ITG iinc);
unsigned long long topodiag_hash_bytes(const void *p,size_t n,
                                       unsigned long long h);
unsigned long long topodiag_hash_d(const double *a,ITG n,
                                   unsigned long long h);
unsigned long long topodiag_hash_i(const ITG *a,ITG n,
                                   unsigned long long h);
unsigned long long topodiag_hash_seed(void);
unsigned long long topodiag_hash_batch(const ITG *elem,ITG n,
                                       ITG *sorted,
                                       unsigned long long h);

double crackcontrol_frame(const double *x1,const double *x2,const double *x3,
                          double *rmat);
double crackcontrol_dir(const double *dl,double beta,double tol,double *m);
double crackcontrol_dissrate(double kn,double tn0,double gc);
ITG crackcontrol_selftest(void);
void crackcontrol_census_zero(crackcontrol_census *s);
/* Eighteen arguments became six.  xstate and v stay parameters on purpose:
   the three call sites pass DIFFERENT arrays - the committed state and
   vold at one, the increment-start xstateini and vini at the other two -
   so they are the caller's choice, not the context's. */
ITG crackcontrol_build(double *c,ITG mode,const trialctx *mdl,
                       const double *xstate,const double *v,
                       crackcontrol_census *s);

#endif /* CCX_FORK_H */
