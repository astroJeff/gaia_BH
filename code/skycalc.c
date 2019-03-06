/* SKY CALCULATOR PROGRAM
   John Thorstensen, Dartmouth College.
   This program computes many quantities frequently needed by the
   observational astronomer.  It is written as a completely
   self-contained program in standard c, so it should be
   very transportable; the only issue I know of that really affects
   portability is the adequacy of the double-precision floating
   point accuracy on the machine.  Experience shows that c compilers
   on various systems have idiosyncracies, though, so be sure
   to check carefully.

   This is intended as an observatory utility program; I assume the
   user is familiar with astronomical coordinates and nomenclature.
   While the code should be very transportable, I also
   assume it will be installed by a conscientious person who
   will run critical tests before it is released at a new site.
   Experience shows that some c compilers generate unforseen errors
   when the code is ported, so the output should be checked meticulously
   against data from other sites.


   After prompting the user for site information, the program enters
   its main command loop,

	- enter RA, dec, proper motion, epoch, date, time,
	     new site parameters, etc. ... (when the code starts,
	     it initializes these things to the current date and time
             and the celestial coordinates of the zenith at the
	     specified site).

	- print an almanac of rise/set and twilight for the sun and
             moon for the specified night

	- compute and display circumstances of observation for the
	   current parameters, including precessed coordinates,
	   airmass, interference from moon or twilight, parallactic
	   angle, etc; the program also gives calendar date in
	   both UT and local, Julian date, and barycentric corrections.

	- compute and display a table of airmasses (etc) at
	   hourly intervals through the night.  This is very useful
	   at the telescope.  Also, if one has a modest number of
	   objects, it may be convenient (using system utilities)
	   to redirect the output and print a hard copy of these
	   tables for ready reference.

	- compute and display galactic and ecliptic coordinates.

	- compute and display rough (of order 0.1 degree, but often
	  much better) positions of the major planets.

    The program is self-contained.  It has been developed primarily on
   UNIX and Linux machines, and should adapt easily to any system with
   a c compiler.

	** BUT CAUTION ... **
   Because many of the routines take a double-precision floating point
   Julian Date as their time argument, one must be sure that the machine
   and compiler carry sufficient mantissa to reach the desired accuracy.
   On most architectures the double-precision floating point julian date
   has an accuracy of order 0.01 seconds of time, which is just adequate.

LEGALITIES:

   I make no guarantee as to the accuracy, reliability, or
   appropriateness of this program, though I have found it to be
   reasonably accurate and quite useful to the working astronomer.

   The program is COPYRIGHT 2000 BY JOHN THORSTENSEN.
   Permission is hereby granted for non-profit scientific or educational use.
   For-profit use (e. g., by astrologers!) must be through negotiated
   license.  The author requests that observatories and astronomy
   departments which install this as a utility notify the author
   by paper mail, just so I know how widely it is used.

   Credits:
    * The julian date and sidereal time routines were
    originally coded in PL/I by  Steve Maker of Dartmouth College.
    They were based on routines in the old American Ephemeris.
    Many of the routines were coded from Jean Meeus' "Astronomical
    Formulae for Calculators", published by Willman-Bell.  This is
    an extraordinarily helpful little book!

    BUGS:

    The program can enter an infinite loop of prompting in response
    to incorrect input.  This is fairly rare in practice, but
    may be bothersome on multi-user systems.  If anyone can document
    this properly, I may be able to fix it.  In the meantime you
    may want to run this only on systems on which you can easily
    stop runaway processes (e.g., by CRTL-C on Linux or UNIX machines.)

    APOLOGIES/DISCLAIMER:
    I am aware that the code here does not always conform to
    the best programming practices.  Not every possible error condition
    is anticipated, and no guarantee is given that this is bug-free.
    Nonetheless, most of this code has been shaken down at several
    hundred sites for several years, and I have never received any
    actual bug reports.  Many users have found this program
    to be useful.

    CHANGES SINCE THE ORIGINAL DISTRIBUTION ....

	The program as listed here is for the most part similar to that
	posted on the IRAF bulletin board in 1990.  Some changes
	include:

	01 In the original code, many functions returned structures, which
	   some c implementations do not like.  These have been eliminated.

	02 The original main() was extremely cumbersome; much of it has
	   been broken into smaller (but still large) functions.

	03 The hourly airmass includes a column for the altitude of the
	   sun, which is printed if it is greater than -18 degrees.

	04 The planets are included (see above).  As part of this, the
	   circumstances calculator issues a warning when one is within
	   three degrees of a major planet.  This warning is now also
	   included in the hourly-airmass table.

	05 The changeover from standard to daylight time has been rationalized.
	   Input times between 2 and 3 AM on the night when DST starts (which
	   are skipped over and  hence don't exist) are now disallowed; input
	   times between 1 and 2 AM on the night when DST ends (which are
	   ambiguous) are interpreted as standard times.  Warnings are printed
	   in both the almanac and calculator mode when one is near to the
	   changeover.

	06 a much more accurate moon calculation has been added; it's used
	   when the moon's coordinates are given explicitly, but not for
	   the rise/set times, which iterate and for which a lower precision
	   is adequate.

	07 It's possible now to set the observatory elevation; in a second
	   revision there are now two separate elevation parameters specified.
	   The elevation above the horizon used only in rise/set calculations
	   and adjusts rise/set times assuming the parameter is the elevation
	   above flat surroundings (e. g., an ocean).  The true elevation above
	   sea level is used (together with an ellipsoidal earth figure) in
	   determining the observatory's geocentric coordinates for use in
	   the topocentric correction of the moon's position and in the
	   calculation of the diurnal rotation part of the barycentric velocity
	   correction.  These refinements are quite small.

	08 The moon's altitude above the horizon is now printed in the
	   hourly airmass calculation; in the header line, its illuminated
	   fraction and angular separation from the object are included,
	   as computed for local midnight.

	09 The helio/barycentric corrections have been revised and improved.
	   The same routines used for planetary positions are used to
	   compute the offset from heliocentric to solar-system
	   barycentric positions and velocities.  The earth's position
	   (and the sun's position as well) have also been improved somewhat.

	10 The printed day and date are always based on the same truncation
	   of the julian date argument, so they should now always agree
	   arbitrarily close to midnight.

	11 A new convention has been adopted by which the default is that the
	   date specified is the evening date for the whole night.  This way,
	   calculating an almanac for the night of July 3/4 and then specifying
	   a time after midnight gives the circumstances for the *morning of
	   July 4*.  Typing 'n' toggles between this interpretation and a
	   literal interpretation of the date.

	12 The planetary proximity warning is now included in the hourly airmass
	   table.

	13 A routine has been added which calculates how far the moon is from
	   the nearest cardinal phase (to coin a phrase) and prints a
	   description.  This information is now included in both the almanac
	   and the calculator mode.

	14 The output formats have been changed slightly; it's hoped this
	   will enhance comprehensibility.

	15 A substantial revision affecting the user interface took place
	   in September of 1993.  A command 'a' has been added to the
	   'calculator' menu, which simply prints the almanac (rise, set,
	   and so on) for the current night.  I'd always found that it was
	   easy to get disoriented using the '=' command -- too much
	   information about the moment, not enough about the time
	   context.  Making the almanac info *conveniently* available
	   in the calculator mode helps your get oriented.

	   When the 'a' almanac is printed, space is saved over the
	   almanac printed on entry, because there does not need
	   to be a banner introducing the calculator mode.  Therefore some
	   extra information is included with the 'a' almanac; this includes
	   the length of the night from sunset to sunrise, the number of
	   hours the sun is below -18 degrees altitude, and the number of hours
	   moon is down after twilight.  In addition, moonrise and moonset
	   are printed in the order in which they occur, and the occasional
	   non-convergence of the rise/set algorithms at high latitude are
	   signalled more forcefully to the user.

	16 I found this 'a' command to be convenient in practice, and never
	   liked the previous structure of having to 'quit' the calculator
	   mode to see almanac information for a different night.  D'Anne
	   Thompson of NOAO also pointed out how hokey this was, especially the
	   use of a negative date to exit. So, I simply removed the outer
	   'almanac' loop and added a 'Q' to the main menu for 'quit'.  The
	   use of upper case -- for this one command only --  should guard
	   against accidental exit.

	17 The menu has been revised to be a little more readable.

	18 More error checking was added in Nov. 1993, especially for numbers.
	   If the user gives an illegal entry (such as a number which isn't
	   legal), the rest of the command line is thrown away (to avoid
	   having scanf simply chomp through it) and the user is prompted
	   as to what to do next.  This seems to have stopped all situations
	   in which execution could run away.  Also, typing repeated carriage
	   returns with nothing on the line -- which a frustrated novice
	   user may do because of the lack of any prompts -- causes a
	   little notice to be printed to draw attention to the help and menu
	   texts.

	19 I found in practice that, although the input parameters and
	   conditions are deducible *in principle* from such things as the
	   'a' and '=' output, it takes too much digging to find them.  So
	   I instituted an 'l' command to 'look' at the current parameter
	   values.  To make room for this I put the 'Cautions and legalities'
	   into the 'w' (inner workings) help text.  This looks as though
	   it will be be very helpful to the user.

	20 The more accurate moon calculation is used for moonrise and
	   moonset; the execution time penalty appears to be trivial.
	   Low precision moon is still used for the summary moon information
	   printed along with the hourly airmass table.

	21 A calculation of the expected portion of the night-sky
	   brightness due to moonlight has been added.  This is based on
	   Krisciunas and Schaefer's analytic fits (PASP, 1991).  Obviously,
	   it's only an estimate which will vary considerably depending on
	   atmospheric conditions.

	22 A very crude estimator of the zenith sky brightness in twilight
	   has been added.

	23 A topocentric correction has been added for the sun, in anticipation
	   of adding eclipse prediction.

	24 The code now checks for eclipses of the sun and moon, by making
	   very direct use of the predicted positions.  If an eclipse is
	   predicted, a notice is printed in print_circumstances; also, a
	   disclaimer is printed for the lunar sky brightness if a lunar
	   eclipse is predicted to be under way.

	25 In the driver of the main calculator loop, a provision has been
	   added for getting characters out of a buffer rather than reading
	   them directly off the command line.  This allows one to type any
	   valid command character (except Q for quit) directly after a number
	   in an argument without generating a complaint from the program
	   (see note 18).  This had been an annoying flaw.

	26 In 1993 December/1994 January, the code was transplanted
	   to a PC and compiled under Borland Turbo C++, with strict
	   ANSI rules.  The code was cut into 9 parts -- 8 subroutine
	   files, the main program, and an incude file containing
	   global variables and function prototypes.

	27 An "extra goodies" feature has been added -- at present it
	   computes geocentric times of a repeating phenomenon as
	   viewed from a site.  This can be used for relatively little-used
           commands to save space on the main menu.

	28 The sun and moon are now included in the "major planets"
	   printout.  This allows one to find their celestial positions
	   even when they are not visible from the selected site.

	29 A MAJOR new feature was added in February 1994, which computes
           the observability of an object at new and full moon over a
           range of dates.  The galactic/ecliptic coordinate converter
           was moved to the extra goodies menu to make room for this.

	30 Inclusion of a season-long timescale means that it's not
           always necessary to specify a date on entry to the program.
           Accordingly, the program immediately starts up in what used
           to be called "calculator" mode -- only the site is prompted
           for.  It is thought that the site will be relevant to nearly
           all users.

	31 Because the user is not led by the hand as much as before, the
           startup messages were completely revised to direct new users
           toward a short `guided tour' designed to show the program's
	   command structure and capabilities very quickly.  Tests on
	   volunteers showed that users instinctively go for anything
	   called the `menu', despite the fact that that's a slow way to
	   learn, so all mention of the menu option is removed from the
	   startup sequence; they'll catch on soon enough.

	32 Code has been added to automatically set the time and
           date to the present on startup.  A menu option 'T' has been
           added to set the time and date to the present plus a settable
           offset.  This should be very useful while observing.

	33 Because Sun machines apparently do not understand ANSI-standard
           function declarations, the code has been revised back to K&R
           style.  It's also been put back in a monolithic block for
           simplicity in distribution.

	34 The startup has been simplified still further, in that the
           coordinates are set automatically to the zenith on startup.
	   An 'xZ' menu item sets to the zenith for the currently specified
           time and date (not necessarily the real time and date.)

	35 Another MAJOR new capability was added in early 1994 --
           the ability to read in a list of objects and set the current
	   coordinates to an object on the list.  The list can be sorted
           in a number of ways using information about the site, date
           and time.

	35 Calculator-like commands were added to the extra goodies menu
           to do precessions and to convert Julian dates to calendar
           dates.  An option to set the date and time to correspond to
           a given julian date was also added.

	36 Another substantial new capability was added Aug 94 -- one can
           toggle open a log file (always named "skyclg") and keep
           a record of the output.  This is done simply by replacing
           most occurrences of "printf" with "oprintf", which mimics
           printf but writes to a log file as well if it is open.
	   This appears to slow down execution somewhat.

	37 12-degree twilight has been added to the almanac.  While the
	   awkward "goto" statements have been retained, the statement
           labels have been revised to make them a little clearer.

	38 The precession routine was generalized to include nutation and
	   aberration, and routines to calculate these effects are now
           included.  Nearly all calls to the new routine leave the
           aberration and nutation out, but 'xa' extra goodies gives
           apparent place up to aberration and nutation.  Tests against the
	   FK5 and "Apparent Places of Fundamental Stars" shows the
           agreement to 0.1 arcsec or smaller.

	39 Atmospheric refraction is computed and reported in the "xa"
           apparent place calculation.  Barometric pressure is guesstimated
           from the observatory elevation above sea level.

	40 An airmass based on a series expansion is printed in the "="
           command, in place of the secant of the zenith distance, provided
           sec z is < 12.   The approximation breaks down larger than this.

	41 The parallactic angle calculation has been simplified and folded
           inside the routine which calculates altitude and azimuth.

	42 A few utility routines (atan_circ, xyz_cel) were cleaned up.

	43 I added a new option for reading the system clock every time
	   'timely' output is asked for, and computing for the updated
	   time & date.  This is toggled on and off with 'xU', for
	   updating.  The date and time can optionally be offset from that
	   read from the system clock, as in the 'T' option.

	44 The sexigesimal coordinate-output routine "put_coords" was
           revised to be more general and a little less ugly.  It's not
           necessarily shorter!  This required gussying up the outputs
           because the new routine doesn't space just the same.
           It's now possible to include a "+" sign explicitly in the
           output when desired.

        45 The calendrical conversion routines have been replaced with
           routines coded from Meeus' "Astronomical Formulae for
           Calculators", to avoid the use of proprietary Numerical
           Recipes code.

        46 The routine used to sort objects by airmass, etc. has been
           replaced with an original coding of the heapsort algorithm
           to avoid the use of proprietary "Numerical Recipes" code.
           All the code *should* now be freely distributable.

        47 The routine which reads strings from input has been rewritten
           so that it can take delimited strings like 18:02:12.33 ...
           delimiter can be any character which is not a numeral,
           a + or - sign, or a decimal point (i.e., anything which can't
           be part of a number.)

	48 Another 'xtra goodies' option has been added to print out the
	   parallax factors (i.e., parallax displacement in ra and dec
           a star would have at exactly 1 pc distance at the nominal
           time and date).  Also prints XYZ for earth and annual aberration
           (displacement of star due to finite speed of light and earth's
           motion about the sun.)

        49 In 2000 October I fixed some infelicities with the
           xv and xf commands, and updated the xd (TDT - UT) command.
           The TeX documentation and the on-line documentation were
	   updated.  In 2002 Jan. finally fixed a bug which displayed the
	   wrong ephemeris in the header in the xv command
	   (computations were OK).

	50 Added the capability to read and write site-parameter files.

	51 replaced "round" with "roundx" to avoid name conflict arising
	   with a math library function.

	52 2003 -- updated the etcorr routine again to reflect Delta T
	   values in 2003 almanac and adjust forward extrapolation to
	   match last tabulated value.

        53 Added code adapted from F. Ochsenbein (CDS Strasbourg) to give the
	   constellation location of the point in question, and report the
	   location in the "=" and "xc" commands.

	54 Added a "parallactic penalty" factor, i.e. the magnitude of
	   tan z times the magnitude of the sine of the parallactic angle.
	   This severity of the parallactic differential refraction effect
	   will be proportional to this, for a slit which defaults to
	   north-south.

	55 Updated the Daylight Saving Time prescription for the U.S. in
	   accordance with the new energy bill.  This applies only from
	   2007 onward.  Also, changed incorrect references to "Daylight
	   Savings" (a mistake noted by F. A. Ringwald).

*/

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <stdarg.h>
#include <string.h>
#include "skycalc.h"


// Define globals
int update_on = 0;
double update_delta = 0.;
FILE *sclogfl = NULL;


main()

{
	struct date_time date,dateback;
	struct coord ttime;
	double jd, jdmid, jdc, jdtmp;
	double jdb, jde, test;  /* jd of begin and end of dst */
	double sid, sss;
	double Dt; /* ephemeris time correction */
	short option, trying;
	short done = 0, optdone = 0, valid_date = 0, nreturns=0, nxret = 0;
	short day, yearday, dst=0, dow;
	int cc, end_in;  /* control character for circumstances loop */
	int cx;    /* control character for extra goodies ... */
        char cxstr[3];
	double objra=0., objdec=0., objepoch=2000.,dectest, ep;
	char constelname[5];   /* 3-letter constellation abbreviation. */
	double curra, curdec, curep, obj_moon;
	double pra[10],pdec[10];
	double alt, az, ha, secz, jdloc;
	double rasun,decsun,ramoon,decmoon,distmoon;
	short enter_ut = 0; /* are times to be entered as UT? */
	short night_date = 1; /* interperet current date as evening or true? */
	char str[80]; /* dummy string */
	char errprompt[80];
	short nch;
	double glong, glat, eclong, eclat, par;
	int status;
	double mura_sec=0.,mura_arcs=0.,mudec=0.,mura,mudc;  /* proper motions */
	double objra_adj, objdec_adj;           /* equinox of std epoch,
			adjusted for proper motion only */
	short hr_span,i;  /* for table of hour angles */
	double tcor, vcor, vcorlsr; /* time and velocity helio corrections */
	char obs_code;
	double minoffset; /* minutes offset from system clock ... */
        double aberra, aberdec, parra, pardec;  /* aberration and
		parallax factors ... */

	/* all the site-specific quantities are here:
		longit     = W longitude in decimal hours
		lat        = N latitude in decimal degrees
		stdz       = standard time zone offset, hours
		use_dst    = 1 for USA DST, 2 for Spanish, negative for south,
				 0 to use standard time year round
		zone_name  = name of time zone, e. g. Eastern
		zabr       = single-character abbreviation of time zone
		site_name  = name of site.
	        elevsea    = elevation above sea level, meters
		elev       = elevation above local horizon, meters
	*/

	/* Kitt peak, MDM observatory, is initialized here as a default.
	   User later gets to choose from a menu of possible sites -
	   they're all hard-coded in the routine 'load_site'. */

	char site_name[45];  /* initialized later with
			   strcpy for portability */
	char zabr = 'M';
	char zone_name[25]; /* this too */
	short use_dst = 0;
	double longit = 7.44111;
	double elevsea = 1925.;  /* for MDM, strictly */
	double elev = 500.; /* well, sorta -- height above horizon */
	double horiz = 0.7174;
	double lat = 31.9533;
	double stdz = 7.;

	/* and off we go.... */

	strcpy(site_name,"Kitt Peak");
	strcpy(zone_name,"Mountain");

	star_tzero = 0.;
	star_period = 0.;
	star_terr = 0.;
	star_perr = 0.;

	oprntf("\nAstronomical calculator program, by John Thorstensen.\n\n");

	load_site(&longit,&lat,&stdz,&use_dst,zone_name,&zabr,
			&elevsea,&elev,&horiz,site_name);
	oprntf("You have selected %s\n",site_name);
        oprntf("(This can be changed later using the 's' [site] command.)\n\n");

#if SYS_CLOCK_OK == 1

        if(get_sys_date(&date,use_dst,enter_ut,night_date,stdz,0.) != 0) {
	  date.y = 2000;  /* have to have a default date.*/
	  date.mo = 1;
	  date.d = 1;
	  date.h = 0.;
	  date.mn = 0.;
	  date.s = 0.;
	  oprntf("SYSTEM CLOCK didn't read. Time & date set arbitrarily to\n");
	  print_all(date_to_jd(date));
	  oprntf("\n");
        }

        else set_zenith(date, use_dst, enter_ut, night_date, stdz, lat,
	  longit, objepoch, &objra, &objdec);

#else
       	  date.y = 2000;  /* have to have a default date.*/
	  date.mo = 1;
	  date.d = 1;
	  date.h = 0.;
	  date.mn = 0.;
	  date.s = 0.;
	  oprntf("SYSTEM CLOCK options turned off, so \n ");
	  oprntf("time and date set arbitrarily to:\n");
	  print_all(date_to_jd(date));
	  oprntf("\n\n");
#endif

        oprntf("\nREADY TO COMPUTE.  Use simple commands to set the date, time, RA\n");
        oprntf("dec, and so on; then use other commands to compute and display\n");
        oprntf("observability information.\n\n");
        oprntf("NEW or RUSTY USERS: type 'f' (and return) for FAST guided tour.\n");

	while((cc = getch()) != 'Q')    switch(cc) {
		case '?':    /* print a menu */
			print_menu();
			nreturns=0;
			break;
                case 'f':    /* print a short tutorial */
                        print_tutorial();
			nreturns=0;
			break;
		case 'r':   /* enter the object's right ascension */
			objra = get_coord();
			nreturns=0;
			break;
		case 'd':   /* enter the object's declination */
			/* filter declination in put for 'date' input! */
			dectest = get_coord();
			if(fabs(dectest) <= 90.) {
				objdec = dectest;
			}
			else {
			    oprntf("REJECTED 'd' INPUT - DECLINATION MUST BE < 90.\n");
			    oprntf("if you want DATE, Use 'y' (yyyy mm dd)\n");
			}
			nreturns=0;
			break;

		case 'C':  /* enter RA and dec together ... */
			 objra = get_coord();
	                        /* filter declination in put for 'date' input! */
                        dectest = get_coord();
                        if(fabs(dectest) <= 90.) {
                                objdec = dectest;
                        }
                        else {
                            oprntf("REJECTED 'd' INPUT - DECLINATION MUST BE < 90.\n");
                            oprntf("if you want DATE, Use 'y' (yyyy mm dd)\n");
                        }
                        nreturns=0;
                        break;

		case 'p':   /* enter the object's proper motions -- */

                	status = get_pm(objdec,&mura_sec,&mudec);
			nreturns=0;
			break;
		case 'e':   /* enter the input epoch */
			getdouble(&objepoch,-10000.,10000.,
				"Give input epoch again...\n"); /* liberal lims*/
                        if(objepoch < -5000.) {
				objepoch = 2000.+
	                            (true_jd(date, use_dst, enter_ut, night_date, stdz)
				     -J2000)/365.25;
				printf("LARGE NEGATIVE EPOCH --- causes input epoch to be set to current!\n");
				printf("set to Julian epoch %9.4f\n",objepoch);
                        }
			break;
		case 't':   /* enter the time ... hours min sec */
			get_time(&date);
			nreturns=0;
#if SYS_CLOCK_OK == 1
			if(update_on == 1) {
				update_on = 0;
				printf("Automatic updating DISABLED.\n");
			}
#endif
			break;
		case 'T':   /* read system clock -- set date & time to that. */
			;
#if SYS_CLOCK_OK == 1
			printf("Set to how many minutes into the future? :");
			scanf("%lf",&minoffset);
			get_sys_date(&date,use_dst,enter_ut,night_date,
				stdz,minoffset);
			if(update_on == 1) {
				update_on = 0;
				printf("Automatic updating DISABLED.\n");
			}
#else
			printf("Sorry -- system clock options are disabled, probably because of an\n");
			printf("incompatibility between your system and the standard time library\n");
			printf("functions used in the program.  Turning them on would require fixing\n");
			printf("the problem in source code and recompiling.\n");
#endif
			nreturns=0;
                        break;
		case 'g':   /* toggle whether times are entered as
				      Greenwich or local */
			if(enter_ut == 1) enter_ut = 0;
			else {
				enter_ut = 1;
				night_date = 0;
			}
			if(enter_ut == 1)
				oprntf("Dates and times entered are now UT.\n");
else oprntf("Dates & times entered are local, dates are literal (not evening).\n");
			oprntf("TIME IS CHANGED to %d %02d %02d, %02d %02d %02.0f",
				date.y,date.mo,date.d,date.h,date.mn,date.s);
			if(enter_ut == 1) oprntf(" UNIVERSAL time.\n");
				else oprntf(" LOCAL time.\n");
			nreturns=0;
			break;
		case 'n': /* toggle whether the current date is to
				  be interpreted as the evening date (for
				  all night) or the true date .... */
			if(enter_ut == 1) {
oprntf("You're entering times as UT, so 'evening date' makes no sense!....\n");
oprntf("No action taken on 'n', first use 'g' first to enable local time input.\n");
			}
			else if(night_date == 1) {
			   night_date = 0;
oprntf("The date in effect will now be interpreted literally, not as evening.\n");
			}
			else {
			   night_date = 1;
oprntf("The date in effect will now be interpreted as the evening date.\n");
			}
			nreturns=0;
			break;
		case 'a':
			if(sclogfl != NULL) fprintf(sclogfl,"\n\n");
					/* space it */
			oprntf("*** Almanac for the currently specified date ***");
			if(night_date != 1) {
oprntf(", but CAUTION!!\nThe 'night date' option is off, so be especially careful\n");
oprntf("this is the correct night.  See 'g' and 'n'...");
			}
			else oprntf(":");
			/* require a trap here to catch bad dates .... */
			if(date.y < 1901 || date.y > 2099)
				oprntf("Bad date!  Trapped! \n");
			else {
			    print_tonight(date,lat,longit,elevsea,elev,horiz,site_name,stdz,
			       zone_name,zabr,use_dst,&jdb,&jde,2);
			    printf("\nType command, 'f' for fast tour, or '?' for menu:");
			}
			nreturns=0;
			break;
		case 'y':   /* enter the date, yyyy mm dd */
			get_date(&date);
			nreturns=0;
#if SYS_CLOCK_OK == 1
			if(update_on == 1) {
				update_on = 0;
				printf("Automatic updating DISABLED.\n");
			}
#endif
			break;

	/* The site parameters must all be changed at the same time; hence
	   user is forced to change them all. */

		case 's':  /* change the site parameters */
			load_site(&longit,&lat,&stdz,&use_dst,
			   zone_name,&zabr,&elevsea,&elev,
			   &horiz,site_name);
			oprntf("New site = %s\n",site_name);
			printf("(Give command, or ? for menu.)\n");
			nreturns=0;
			break;
		case '=':  /* PRINT CIRCUMSTANCES for current params */
#if SYS_CLOCK_OK == 1
			if(update_on == 1) {
        			if(get_sys_date(&date,use_dst,enter_ut,
					night_date,stdz,update_delta) != 0)
				       printf("Can't get system date! \n");
			}
#endif
			if(sclogfl != NULL) fprintf(sclogfl,"\n\n*** Instantaneous Circumstances ***\n");
			if(setup_time_place(date,longit,lat,stdz,
			    use_dst,zone_name,zabr, site_name,enter_ut,
			    night_date,&jd,&jdloc,&jdb,&jde,&sid,
			    &curep) < 0) break;
			print_circumstances(objra,objdec,objepoch,jd,
			    curep,mura_arcs,mura_sec,mudec,
				   sid,lat,elevsea,horiz);
			nreturns=0;
			break;
		case 'm':  /* print positions of major planets */
#if SYS_CLOCK_OK == 1
			if(update_on == 1) {
        			if(get_sys_date(&date,use_dst,enter_ut,
					night_date,stdz,update_delta) != 0)
				       printf("Can't get system date! \n");
			}
#endif
			if(setup_time_place(date,longit,lat,stdz,
				use_dst,zone_name,zabr, site_name,enter_ut,night_date,
				&jd,&jdloc,&jdb,&jde,&sid,&curep) < 0)
				   break;
			comp_el(jd);
			pposns(jd,lat,sid,1,pra,pdec);
			nreturns=0;
			break;
		case 'h':  /* print an hourly airmass table */
#if SYS_CLOCK_OK == 1
			if(update_on == 1) {
        			if(get_sys_date(&date,use_dst,enter_ut,
					night_date,stdz,update_delta) != 0)
				       printf("Can't get system date! \n");
			}
#endif
			hourly_airmass(date,stdz,lat,longit,horiz,
			   use_dst,objra,objdec,objepoch, mura_sec,
			   mura_arcs,mudec);
			nreturns=0;
			break;
		case 'o':
			if(sclogfl != NULL) fprintf(sclogfl,"\n");
		        obs_season(objra,objdec,objepoch,
			     lat,longit);
			nreturns=0;
                        break;
		case 'c':  /* print galactic and ecliptic coordinates */
			galact(objra,objdec,objepoch,&glat,&glong);
			oprntf("Galactic: l = %5.2f, b = %5.2f\n",
				glat,glong);
			eclipt(objra,objdec,objepoch,date_to_jd(date),
				&curep,&eclong,&eclat);
			oprntf("Ecliptic (equinox %7.2f): long = %5.2f, lat = %5.2f\n",
				curep, eclong, eclat);
                        radec_to_constel(objra, objdec, objepoch, constelname);
			oprntf("Located in contellation %s.\n",constelname);

			nreturns=0;
			break;
		case 'l':
			print_params(date,enter_ut,night_date,
				stdz,lat,longit,site_name,elevsea,elev,use_dst,
				objra,objdec,objepoch,mura_sec,mura_arcs,mudec);
			nreturns=0;
			break;
		case 'i':  /* print a short tutorial */
			print_examples();
			nreturns=0;
			break;
		case 'w':  /* print information about algorithms, acc. */
			print_accuracy();
			nreturns=0;
			break;
		case 'x':
			nxret = 0;
    /*			printf("(Give xtra goodies subcommand, ? for menu)\n");
			scanf("%s",cxstr);    */
		        while(isspace(cx = getch()) != 0) {
			    nxret++;
			    if(nxret == 3) {
				printf("Give an extra goodies command, or ? for menu!\n");
				nxret = 0;
			    }
			}
            /*          cx = cxstr[0];   */
			switch(cx)  {
                           case '?':
			oprntf("Extra goodies commands are:\n");
			oprntf("  x? ... print extra goodies menu.\n");
			oprntf("  xy ... print day of year.\n");
			oprntf("  xc ... give galactic and ecliptic coords & lsr corr'n.\n");
			oprntf("  xa ... print apparent place (nutation and aberration corrected)\n");
			oprntf("  xd ... give rough value of delta T = TDT - UT.\n");
                        oprntf("  xv ... list geocentric times of repeating phenom (Variable star)\n");
			oprntf("  xf ... give phase of repeating phenom.\n");
			oprntf("  xp ... parallax factors and aberration.\n");
                        oprntf("  xb ... precess a bunch of coords the same way.\n");
			oprntf("  xj ... calculate calendar dates given julian dates.\n");
			oprntf("  xJ ... *set* date and time values from Julian date.\n");
			oprntf("  xZ ... *set* RA and dec to Zenith\n");
#if SYS_CLOCK_OK == 1
			oprntf("  xU ... *toggle* automatic system-clock update.\n");
#endif

#if LOG_FILES_OK == 1
			oprntf("LOG-FILE COMMAND:\n");
			oprntf("  xL ... toggles log file open or closed\n");
#endif
			oprntf("COMMANDS FOR FILES OF OBJECTS:\n");
			oprntf("  xR ... read objects from a file, format: name h m s d m s epoch\n");
			oprntf("  xl ... type out (part of) object list.\n");
			oprntf("  xN ... find object by name, set to its coords\n");
			oprntf("  xS ... sort and select object by a rank, set to coords.\n");
			oprntf("(Note that capital letters affect more than one quantity, eg. RA and dec)");
			   break;

                           case 'v':
			      if(sclogfl != NULL)
				fprintf(sclogfl,"\n\n  *** Ephemeris predictions ***\n\n");
			      ephemgen(objra,objdec,objepoch,lat,longit);
			   break;

			   case 'f':
#if SYS_CLOCK_OK == 1
			      if(update_on == 1) {
        			    if(get_sys_date(&date,use_dst,enter_ut,
					night_date,stdz,update_delta) != 0)
				       printf("Can't get system date! \n");
			      }
#endif
/*                              print1phase(date, use_dst, enter_ut,
				night_date, stdz, lat,
	  			longit, objepoch, objra, objdec);
*/
				phaselisting(use_dst, enter_ut, night_date, 						stdz, lat, longit,
					objepoch, objra, objdec);
			   break;

			   case 'b':
                              mass_precess();
			   break;

			   case 'a':
#if SYS_CLOCK_OK == 1
			      if(update_on == 1) {
        			    if(get_sys_date(&date,use_dst,enter_ut,
					night_date,stdz,update_delta) != 0)
				       printf("Can't get system date! \n");
			      }
#endif
			      jd = true_jd(date, use_dst, enter_ut, night_date, stdz);
			      print_apparent(objra,objdec,objepoch,mura_sec,mudec,jd,
				lat,longit,elevsea);
			   break;

			   case 'p':
                                 jd = true_jd(date, use_dst, enter_ut,
					night_date, stdz);
				 parellipse(jd, objra, objdec, objepoch, lat,
					 longit, &parra,
                			 &pardec, &aberra, &aberdec);

			   break;

			   case 'c':
			      oprntf("Equatorial: RA = ");
			      put_coords(objra,3,0);
			      oprntf(", dec = ");
			      put_coords(objdec,2,1);
			      oprntf(" (epoch %6.1f)\n",
                                  objepoch);
			      galact(objra,objdec,objepoch,&glat,&glong);
			      oprntf("  Galactic: l = %5.2f, b = %5.2f\n",
				glat,glong);
			      eclipt(objra,objdec,objepoch,date_to_jd(date),
		           		&curep,&eclong,&eclat);
			      oprntf("Ecliptic (equinox %7.2f): long = %5.2f, lat = %5.2f\n",
				curep, eclong, eclat);
			      lsrcor(objra,objdec,objepoch,&vcorlsr);
			      oprntf("Rough correction from helio to lsr %5.1f km/s\n",vcorlsr);
			   break;
			case 'y': /* print day of year ... */
#if SYS_CLOCK_OK == 1
				if(update_on == 1) {
       					if(get_sys_date(&date,use_dst,enter_ut,
						night_date,stdz,update_delta) != 0)
			       		printf("Can't get system date! \n");
				}
#endif

				jdtmp = true_jd(date, use_dst, enter_ut, night_date, stdz);

				printf("\nUT date & time: ");
				print_all(jdtmp);
				printf("\n");
				printf("UT day of year:%10.5f\n",day_of_year(jdtmp));
				break;

			   case 'd':
			      jd = true_jd(date, use_dst, enter_ut, night_date, stdz);
       	                      Dt = etcorr(jd);
			      oprntf("Delta t = TDT - UT = %5.1f seconds\n",Dt);
			      oprntf("JD %f (UT) --> ",jd);
			      jd = jd + Dt / SEC_IN_DAY;
			      oprntf(" %f (TDT)\n",jd);
			      if(date.y > 2001)
				oprntf("(Value is an extrapolated guess ... only computable after the fact.)\n");
			      else oprntf("+- 0.5 sec, based on 5-year linear interpolations.\n");
		              if(date.y < 1983) oprntf("Before 1983, TDT was preceded by ephemeris time (ET).\n");
			   break;
#if SYS_CLOCK_OK == 1
			   case 'U':
				if(update_on == 0) {
					update_on = 1;
					printf("Automatic update toggled ON.\n");
					printf("Offset (minutes) into future?:");
					scanf("%lf",&update_delta);
				}
			        else {
					update_on = 0;
					printf("Automatic update toggled OFF.\n");
				}
		           break;
#endif

			   case 'j':
			      jdc = 1;
			      oprntf("jd to calendar conversion. \n");
			      while(jdc > 0.) {
				  printf("Give jd to convert, negative to exit:");
				  getdouble(&jdc,-1000000.,3000000.,
	  			    "Give JD to convert, negative to exit");
			   	  oprntf("%f -- > ",jdc);
			          print_all(jdc);
			          oprntf("\n");
			      }
			       oprntf("(Value of date in main program is unaffected.)\n");
	                   break;

			   case 'J':
 			      oprntf("Sets date and time from an input JD.\n");
			      printf("Give JD to set to, negative value for no action: ");
				  getdouble(&jdc,-1000000.,LASTJD,
	  			    "Give jd to set to, negative for no action");
			      set_to_jd(&date, use_dst, enter_ut,
					night_date, stdz, jdc,1);
#if SYS_CLOCK_OK == 1
				if(update_on == 1) {
				       update_on = 0;
				       printf("Automatic updating DISABLED.\n");
				}
#endif
	                   break;
#if LOG_FILES_OK == 1
			   case 'L':
				if(sclogfl == NULL) {
				    trying = 1;
				    while(sclogfl == NULL && trying == 1) {
				    	printf("Give filename for log file, type NONE to cancel:");
				    	scanf("%s",str);
					if(strcmp(str,"NONE") == 0)
						trying = 0;
					else {
				    		sclogfl = fopen(str,"a");
					}
				    }
				    if(sclogfl != NULL)
					printf("log file %s is OPEN in append mode.\n",str);
				    else printf("LOG FILE NOT OPENED.\n");
				}
			        else {
				    fclose(sclogfl);
				    sclogfl = NULL;  /* reset it explicitly */
				    printf("Log file has been CLOSED.\n");
				}
				break;
#endif
			   case 'R':
				read_obj_list();
				break;
			   case 'l':
				type_list(date,use_dst,enter_ut,night_date,
				    stdz,lat,longit);
     				break;
			   case 'N':
				find_by_name(&objra,&objdec,objepoch,date,
				  use_dst,enter_ut,night_date,stdz,lat,longit);
				break;
			   case 'S':
				find_nearest(&objra,&objdec,objepoch,date,
                                  use_dst,enter_ut,night_date,stdz,lat,longit);
				break;
			   case 'Z':
                                set_zenith(date,use_dst,enter_ut,night_date,
 					stdz,lat,longit,objepoch,&objra,
					&objdec);
				break;
			   case 't':  /* test */

	test = true_jd(date, use_dst, enter_ut, night_date, stdz);
	printf(".... true_jd gives --> %f\n",test);
	break;
			   case ' ':  ;
			   break;
			   case '\n': ;
    			   break;
			   default: oprntf("Unrecognized character %c ... no action.\n",cx);
                           break;
                        }      /* end of 'extra goodies' menu. */
			nreturns=0;
			printf("\n(eXtra goodies doesn't loop.)\n");
			printf("Back in main commands, 'Q' quits, '?' menu, 'f' fast tour.\n");
			break;
		case '\n': /* ignore carriage returns */
			nreturns++;  /* but guide the user if they keep
			       hitting returns .... */
			if(nreturns == 3) {
			  printf("You're repeating carriage returns. There are no prompts.\n");
			  printf("Type 'f' for fast tour, 'i' for instructions, ? for a menu.\n");
			  nreturns = 0;
			}
			break;
		case ' ':  /* ignore blank spaces */
			break;
		case 'q':  /* prompt if user's trying to quit */
			printf("Type an UPPER CASE Q to quit.\n");
			break;
		default:   /* complain if unrecognizable */
			printf("Unknown command, %c\n",cc);
	}       /* closing switch loop */
	BLUNDER:; /* DUMMY STATEMENT */
	oprntf("Suggestions or comments --> john.thorstensen@dartmouth.edu\n");
	oprntf("Goodbye.\n");
}
