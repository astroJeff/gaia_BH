/* a couple of the system-dependent magic numbers are defined here */

#define SYS_CLOCK_OK 1    /* 1 means ANSI-standard time libraries do work,
   2 means they don't.  This is used by compiler switches in file 5 and
   the main program.  */

#define LOG_FILES_OK 1  /* 1 means that log files are enabled.
			Any other value means they're not.  */

#define MAX_OBJECTS 2000
#define MINSHORT -32767   /* min, max short integers and double precision */
#define MAXSHORT 32767
#define MAXDOUBLE 1.0e38
#define MINDOUBLE -1.0e38
#define BUFSIZE 150

#define XFORM_FROMSTD  1  /* defined quantities for apparent place transforms .. */
#define XFORM_TOSTDEP  -1
#define XFORM_JUSTPRE  1
#define XFORM_DOAPPAR  0
#define PRINT_TO_SCREEN 0

/* some (not all) physical, mathematical, and astronomical constants
   used are defined here. */

#define  PI                3.14159265358979
#define  TWOPI             6.28318530717959
#define  PI_OVER_2         1.57079632679490  /* From Abramowitz & Stegun */
#define  ARCSEC_IN_RADIAN  206264.8062471
#define  DEG_IN_RADIAN     57.2957795130823
#define  HRS_IN_RADIAN     3.819718634205
#define  KMS_AUDAY         1731.45683633   /* km per sec in 1 AU/day */
#define  SPEED_OF_LIGHT    299792.458      /* in km per sec ... exact. */
#define  SS_MASS           1.00134198      /* solar system mass in solar units */
#define  J2000             2451545.        /* Julian date at standard epoch */
#define  SEC_IN_DAY        86400.
#define  FLATTEN           0.003352813   /* flattening of earth, 1/298.257 */
#define  EQUAT_RAD         6378137.    /* equatorial radius of earth, meters */
#define  ASTRO_UNIT        1.4959787066e11 /* 1 AU in meters */
#define  RSUN              6.96000e8  /* IAU 1976 recom. solar radius, meters */
#define  RMOON             1.738e6    /* IAU 1976 recom. lunar radius, meters */
#define  PLANET_TOL        3.          /* flag if nearer than 3 degrees
						to a major planet ... */
#define  KZEN              0.172       /* zenith extinction, mag, for use
				     in lunar sky brightness calculations. */
#define FIRSTJD            2415387.  /* 1901 Jan 1 -- calendrical limit */
#define LASTJD             2488070.  /* 2099 Dec 31 */

/* MAGIC NUMBERS which might depend on how accurately double-
   precision floating point is handled on your machine ... */

#define JDRESOLUTION 3.5e-8   /* 3 milliseconds */

#define  EARTH_DIFF        0.05            /* used in numerical
   differentiation to find earth velocity -- this value gives
   about 8 digits of numerical accuracy on the VAX, but is
   about 3 orders of magnitude larger than the value where roundoff
   errors become apparent. */

#define  MIDN_TOL          0.00001         /* this is no longer
   used -- it was formerly
   how close (in days) a julian date has to be to midnight
   before a warning flag is printed for the reader.  VAX
   double precision renders a Julian date considerably
   more accurately than this.  The day and date are now based
   on the same rounding of the julian date, so they should
   always agree. */

#define ALT_3  19.47  /* 19.47 degrees altitude => sec z = 3 */
#define ALT_2  30.
#define ALT_15 41.81
#define SID_RATE 1.0027379093  /* sidereal / solar rate */

/*  FUNCTION PROTOTYPES and type definitions ....
    These are used in breaking the code into function libraries.
    They work properly on a strictly ANSI compiler, so they
    apparently comply with the ANSI standard format.  */

/* these global variables determine whether the program updates the
time and date whenever asked for output. */

#define ITEMS(a)  (sizeof(a)/sizeof(a[0]))

struct coord
   {
     short sign;  /* carry sign explicitly since -0 not neg. */
     double hh;
     double mm;
     double ss;
   };

struct date_time
   {
	short y;
	short mo;
	short d;
	short h;
	short mn;
	float s;
   };

/* elements of planetary orbits */
struct elements {
   	char name[9];
   	double incl;
   	double Omega;
   	double omega;
   	double a;
   	double daily;
   	double ecc;
   	double L_0;
   	double mass;
   };

struct objct {
	char name[20];
	double ra;
	double dec;
	double mura;
	double mudec;
	float ep;
	float xtra;  /* mag, whatever */
};



void oprntf(char *fmt, ...);
char getch();
void ungetch(int);
int legal_num_part(char);
int legal_int_part(char);
int legal_command_char(char);
int parsedouble(char *, double *);
int getdouble(double *,double,double,char*);
int parseshort(char *, short *);
int getshort(short *, short, short, char *);
double bab_to_dec(struct coord bab);
void dec_to_bab (double, struct coord *bab);
short get_line(char *s);
int is_delimited(char *instrng);
double conv_delimited_coord(char *instrng);
double get_coord();
double roundx(double, int);
void round_coord(struct coord *incoord, struct coord *outcoord, int prec);
void put_hrs(double, short, int, int, int);
void put_coords(double, int, int);
void chomp(char *s);
void load_site(	double *longit, double *lat, double *stdz, short *use_dst,
  char *zone_name, char *zabr, double *elevsea, double *elev, double *horiz,
	char *site_name);
double atan_circ(double, double);
void min_max_alt(double, double, double *, double *);
double altit(double, double, double, double *, double *);
double secant_z(double);
double true_airmass(double);
double ha_alt(double, double, double);
double subtend(double, double, double, double);
int get_pm(double, double *, double *);
int get_date(struct date_time *date);
int get_time(struct date_time *date);
double date_to_jd(struct date_time date);
void caldat(double jdin, struct date_time *date, short *dow);
short day_of_week(double);
double day_of_year(double);
void print_day(short);
void print_all(double);
void print_current(struct date_time date, short night_date, short enter_ut);
void print_calendar(double, short *);
void print_time(double, short);
double frac_part(double);
double lst(double, double);
double adj_time(double x);
void lpmoon(double, double, double, double *, double *, double *);
void lpsun(double, double *, double *);
void eclrot(double, double *, double *, double *);
double circulo(double x);
void geocent(double, double, double, double *, double *, double *);
double etcorr(double jd);
void accumoon(double, double, double, double, double *, double *, double *, double *, double *, double *);
void flmoon(int, int, double *);
float lun_age(double, int *);
void print_phase(double jd);
double lunskybright(double, double, double, double, double, double);
void accusun(double, double, double, double *, double *, double *, double *, double *, double *, double *, double *);
double jd_moon_alt(double, double, double, double, double);
double jd_sun_alt(double, double, double, double);
float ztwilight(double alt);
void find_dst_bounds(int yr, double stdz, int use_dst,
			double *jdb, double *jde);
double zone(short use_dst, double stdz, double jd, double jdb, double jde);
double true_jd(struct date_time date, short use_dst, short enter_ut, short night_date, double stdz);
void print_tz(double jd,short use,double jdb,double jde,char zabr);
void xyz_cel(double x,double y,double z,double* ra,double* dec);
void aberrate(double epoch, double vec[], int from_std);
void nutation_params(double date_epoch, double *del_psi, double *del_ep);
void cooxform(double rin, double din, double std_epoch,
  double date_epoch, double *rout, double *dout, int just_precess, int from_std);
double near_hor_refr(double app_alt, double pressure);
double refract_size(double alt, double elev);
void refract_corr(double *ha, double *dec, double lat, double elev, double *size, int sense);
void mass_precess();
void print_apparent(double, double, double, double, double, double, double, double, double);
void radec_to_constel(double ra, double dec, double epoch, char *constname);
void galact(double, double, double, double *, double *);
void eclipt(double, double, double, double, double *, double *, double *);
void comp_el(double jd);
void planetxyz(int p, double jd, double* x, double *y, double *z);
void planetvel(int p, double jd, double *vx, double *vy, double *vz);
void xyz2000(double, double, double, double);
void xyz2000xf(double, double *, double *, double *);
void earthview(double *x, double *y, double *z, int i, double *ra, double *dec);
void pposns(double jd,double lat, double sid,short print_option,double *planra, double *plandec);
void barycor(double jd,double *x,double *y,double *z,double *xdot,double *ydot,double *zdot);
void helcor(double jd,double ra,double dec,double ha,double lat,double elevsea,double *tcor,double *vcor);
void lsrcor(double ra,double dec,double epoch,double *vcor);
void parellipse(double jd, double ra, double dec, double epoch, double lat,
    double longit, double *dra,
		double *ddec, double *aberra, double *aberdec);
float overlap(double r1, double r2, double sepn);
void solecl(double sun_moon,double distmoon,double distsun);
short lunecl(double, double, double, double, double, double);
void planet_alert(double, double, double, double);
short setup_time_place(struct date_time date, double longit, double lat, double stdz,
  short use_dst, char *zone_name, char zabr, char *site_name, short enter_ut,
  short night_date, double *jdut, double *jdlocal, double *jdb, double* jde,
  double *sid, double *curepoch);
void print_tonight(struct date_time date, double lat, double longit, double elevsea,
  double elev, double horiz,char *site_name, double stdz, char *zone_name, char zabr,
  short use_dst, double *jdb, double *jde, short short_long);
void print_circumstances(double, double, double, double, double, double, double,
  double, double, double, double, double);
void hourly_airmass(struct date_time date, double stdz, double lat, double longit,
  double horiz, short use_dst, double objra, double objdec, double objepoch,
  double mura_sec, double mura_arcs, double mudec);
void print_params(struct date_time date,short enter_ut, short night_date,double stdz,
  double lat, double longit, char *site_name, double elevsea, double elev, short use_dst,
  double objra, double objdec, double objepoch, double mura_sec, double mura_arcs, double mudec);
void print_menu();
void print_tutorial();
void print_examples();
void print_accuracy();
void print_legalities();
double hrs_up(double jdup, double jddown, double jdeve, double jdmorn);
void print_air(double secz, short prec);
void print_ha_air(double ha, double secz, short prec1, short prec2);
void obs_season(double, double, double, double, double);
int get_sys_date(struct date_time *date, short use_dst, short enter_ut, short night_date,
  double stdz, double toffset);
void inheap(float *a, int *ind, int *n);
void setheap(float *a, int *ind, int n);
void swapheap(float *a, int *ind, int n, int inputind, int *outind);
void outheap(float *a, int *ind, int *n, int inputind, int *outputind);
void indexx(int n, float *a, int *ind);
int read_obj_list();
int find_by_name(double *ra, double *dec, double epoch, struct date_time date,
  short use_dst, short enter_ut, short night_date, double stdz, double lat, double longit);
void type_list(struct date_time date, short use_dst, short enter_ut, short night_date,
  double stdz, double lat, double longit);
int find_nearest(double *ra, double *dec, double epoch, struct date_time date, short use_dst,
  short enter_ut, short night_date, double stdz, double lat, double longit);
void set_zenith(struct date_time date, short use_dst, short enter_ut, short night_date,
  double stdz, double lat, double longit, double epoch, double *ra, double *dec);
void print1phase(struct date_time date, short use_dst, short enter_ut, short night_date,
  double stdz, double lat, double longit, double epoch, double ra, double dec);
int set_to_jd(struct date_time *date, short use_dst, short enter_ut, short night_date,
  double stdz, double jd, int verbose);
void phaselisting(short use_dst, short enter_ut, short night_date, double stdz,
  double lat, double longit, double epoch, double ra, double dec);


#if SYS_CLOCK_OK == 1
extern int update_on;
extern double update_delta;
// extern int update_on = 0;
// extern double update_delta = 0.;
#endif


extern FILE *sclogfl;
// FILE *sclogfl = NULL;

extern double star_tzero, star_terr, star_period, star_perr;  /* for ephemeris calculations ... global */
