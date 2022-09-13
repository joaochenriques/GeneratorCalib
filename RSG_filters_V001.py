
#~#############################################################################
# Version: 2010/10/18 17h11
#~#############################################################################
#  LOW-PASS filter
#  [ 15, 13, 11,  9,  7,  5 ]
#  [ 25, 23, 21, 19, 17, 15 ]
#  [ 35, 33, 31, 29, 27, 25 ]
#  [ 39, 37, 35, 33, 31, 29 ]
#  [ 45, 43, 41, 39, 37, 35 ]
#  [ 49, 47, 45, 43, 41, 39 ]
#  [ 61, 59, 55, 53, 51, 47 ]
#  [ 69, 65, 63, 59, 57, 53 ]
#  [ 77, 73, 71, 67, 63, 59 ]
#  [101, 97, 91, 87, 83, 77 ]
#~#############################################################################
#  HIGH-PASS filter
#  [ 177, 169, 161, 153, 145, 135 ]
#  [ 237, 225, 215, 203, 193, 181 ]
#  [ 295, 281, 267, 253, 239, 227 ]
#  [ 355, 339, 321, 305, 289, 273 ]
#  [ 413, 393, 375, 355, 335, 317 ]
#  [ 433, 413, 395, 375, 355, 337 ]
#  [ 453, 443, 425, 395, 385, 367 ]  -> JH
#  [ 463, 453, 435, 405, 395, 377 ]  -> JH
#  [ 473, 451, 429, 407, 385, 363 ]
#  [ 531, 507, 481, 457, 431, 407 ]
#  [ 591, 563, 535, 507, 481, 453 ]
#  [ 581, 553, 527, 499, 473, 445 ]
#  [ 727, 693, 659, 625, 591, 557 ]
#  [ 871, 831, 789, 749, 707, 667 ]
#  [1017, 969, 921, 873, 825, 779 ]
#  [1161,1107,1051, 997, 943, 889 ]
#  [1307,1245,1185,1123,1061, 999 ]
#  [1451,1383,1315,1247,1179,1109 ]
#~#############################################################################
from pylab import plot, legend, show, hold, figure, subplot, semilogy, \
                  xlabel, ylabel, grid, ion, ioff, draw, semilogy, pi, loglog,\
                  title, axis, yticks, close, randn, gcf, average, unwrap, \
                  angle, real, imag
from numpy import concatenate, zeros, ones, array, polyval, arange, log10,\
                  dot, linalg, polyfit, poly1d, linspace, polyval, fliplr

from scipy.signal import freqz

#~#############################################################################
def filter_multiply( filter, signal ):
  from scipy import weave, zeros
  from scipy.weave import converters
  from sys import platform

  half_window = int( filter.size / 2 )

  if filter.size % 2 == 0:
    raise "ERROR@filter_multiply: filter size is even."

  out_signal = zeros( signal.size - 2*half_window )

  if platform.find( 'linux' ) >= 0:

    support = """
    #include <iostream>
    """

    code = """
    int ej = filter.extent(0);
    int hw = ej / 2;
    int ei = signal.extent(0) - 2*hw;

    for( int i = 0; i < ei; i++ )
    {
      for( int j = 0; j < ej; j++ )
        out_signal(i) += signal(i+j) * filter(j);
    }
    """

    weave.inline( code, [ 'signal', 'filter', 'out_signal' ],
                  type_converters = converters.blitz,
                  support_code = support )
  else:
    ej = filter.size;
    hw = ej / 2;
    ei = signal.size - 2*hw;

    for i in range( 0, ei ):
      out_signal[i] = dot( signal[i:i+ej], filter )

  return out_signal

#~#############################################################################
def filters_plot( bhigh, f_sample, numfig, clr='r' ):

  filter_name = []
  first = True
  figure( numfig )

  ax1 = subplot(111)
  ylabel( 'Amplitude' )
  xlabel( 'Frequency (Hz)' )
  grid()

  freqs = linspace( 0.0, 1 / 10.0, 50000. )
  nplot = len( freqs )

  wf = None
  hf = ones( nplot )

  for b in bhigh:
    n = b.size
    filter_name.append( n )
    a = zeros(n)
    a[ n / 2 ] = 1.0
    w, h = freqz( b, a, freqs )
    hf *= abs(h)

    ah  = 20 * log10( abs(h) )

    if first:
      # freq_filter = freq_normalized * freq_sample / ( 2 * pi )
      wf  = w.copy()
      wf *= f_sample / ( 2.0 * pi )

    plot( wf, ah, 'b' if first else 'c' )
    first = False

  angles = unwrap( angle( hf ) )
  hf = 20 * log10( hf )

  plot( wf, hf, clr, label = 'filter' )

  ax2 = ax1.twinx()
  plot( wf, angles, 'm')

  title( ( 'Filter %s @ f_s = %.1f Hz' % ( str( filter_name ), f_sample ) ) )
  legend(loc = 'upper right')
  ylabel('Phase angle (radians)', color='m')

  return

#~#############################################################################
def filters_sym_6th( pnts_seq ):

  b_lst = []
  for Npts in pnts_seq:
    b = zeros( Npts )

    float_Npts = float( Npts )
    hpts = Npts / 2

    if Npts % 2 == 0:
      print( ' Number of points should be an odd ', Npts )
      exit(1)

    CONST6_0 = (   1225.E0*float_Npts**6 -  57575.E0*float_Npts**4 +  \
                 605395.E0*float_Npts**2 - 952245.E0 ) / 64.E0
    CONST6_2 = ( -11025.E0*float_Npts**4 + 330750.E0*float_Npts**2
               -1507485.E0 ) / 16.E0
    CONST6_4 = (  24255.E0*float_Npts**2 - 347655.E0 ) / 4.E0
    CONST6_6 = -15015.E0
    Denomi_6 =  4.E0 * float_Npts * ( float_Npts**6 - 56.E0*float_Npts**4
            + 784.E0*float_Npts**2 - 2304.E0 )

    mid_Coef = array( [ CONST6_0 / Denomi_6 ] )
    pl = array( [CONST6_6, 0.0, CONST6_4, 0.0, CONST6_2, 0.0, CONST6_0] )

    float_I = arange( 1.0, hpts + 0.1, 1.0 )

    Coefs = polyval( pl, float_I ) / Denomi_6


    b_lst.append( concatenate( ( Coefs[::-1], mid_Coef, Coefs ) ) )

  return b_lst

#~#############################################################################
def LS_fitting( time, signal, degree ):
  return poly1d( polyfit( time, signal, degree ) )

#~#############################################################################
def sym_filters_avg_ends( b_lst, y ):
  # average ends
  for b in b_lst:
    half_window = int( len(b) / 2 )

    # average seems to work well -> good first approach
    firstvals = ones( half_window ) * average( y[ 1: half_window+1] )
    lastvals  = ones( half_window ) * average( y[-half_window-1:-1] )

    # mirrors -> does not work very well...
    #~ firstvals = y[ 0] - abs( y[ 1: half_window+1][::-1] - y[ 0] )
    #~ lastvals  = y[-1] - abs( y[-half_window-1:-1][::-1] - y[-1] )

    y = concatenate( ( firstvals, y, lastvals ) )
    y = filter_multiply( b, y )

  return y

#~#############################################################################
def sym_filters_czo_ends( b_lst, y ):
  # continuous zero orders ends
  window = 4*len( b_lst[0] )

  half_window = 0
  for b in b_lst:
    half_window += int( len(b) / 2 )

  firstvals = ones( half_window ) * average( y[ :window] )
  lastvals  = ones( half_window ) * average( y[-window:] )
  y = concatenate( ( firstvals, y, lastvals ) )

  for b in b_lst:
    y = filter_multiply( b, y )

  return y

#~#############################################################################
def sym_filters_cls_ends( b_lst, y, ends_degree = 6, mult = 1 ):
  # continuous least-square ends
  for b in b_lst:
    window      = mult*len(b)
    half_window = int( len(b) / 2 )

    time        = arange( 0, window+0.5, 1.0 )

    left_time   = time[ :half_window ]
    left_y      = y[ : window + 1 ]
    poly_left   = LS_fitting( time, left_y,  ends_degree )
    left_vals   = poly_left ( left_time )

    right_time  = time[ -half_window: ]
    right_y     = y[ -window-1: ]
    poly_right  = LS_fitting( time, right_y, ends_degree )
    right_vals  = poly_right( right_time )

    y = filter_multiply( b, y )
    y = concatenate( ( left_vals, y, right_vals ) )

  return y

#~#############################################################################
def filters_SG2( pnts_seq, degree ):
  b_lst = []
  for Npts in pnts_seq:

    M = Npts / 2
    b = zeros( Npts )
    b[M+1] = 1
    q = linspace(-M,M,Npts)
    #~ print( len(b), len(q) )

    a = polyfit(q, b,degree);

    h = polyval(a,q)

    if Npts % 2 == 0:
      print( ' Number of points should be an odd ', Npts )
      exit(1)

    b_lst.append( h )

  return b_lst

#~#############################################################################
if __name__=='__main__':

  #~ f_sample = 62500.0
  #~ high_filter_seq = array( [ 177, 169, 161, 153, 145, 135 ] )
  #~ high_filter_seq2 = array( [1019, 989, 959, 929, 899, 869 ] )
  #~ low_filter_seq =  array( [2903, 2767, 2631, 2495, 2359, 2219] )

  f_sample = 62500.0
  high_filter_seq = array( [1063, 1015,  967,  919,  871,  811] )
  #~ high_filter_seq = array( [ 727, 693, 659, 625, 591, 557 ] )
  #~ high_filter_seq = array( [6115, 5885, 5705, 5525, 5369] )
  low_filter_seq =  array( [17419, 16603, 15787, 14971, 14155, 13315] )

  #~ high_filter_seq  *= 6
  #~ high_filter_seq  += 1

  #~ high_filter_seq2 *= 6
  #~ high_filter_seq2 += 1

  #~ low_filter_seq   *= 6
  #~ low_filter_seq   += 1

  high_b = filters_sym_6th( high_filter_seq )
  #~ high_b2 = filters_SG2( high_filter_seq2, 86 )

  #~ bw = filters_plot( high_b, f_sample, 30 )
  bw = filters_plot( high_b, f_sample, 30, 'k' )

  low_b = filters_sym_6th( low_filter_seq )
  bw = filters_plot( low_b, f_sample, 31 )

  show()

#~EOF##########################################################################
