/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package PedigreeSim;

import java.util.ArrayList;
import java.util.Random;

/**
 * Tools provides some generic routines
 * One instance of Tools needs to be constructed, because Tools also provides
 * a Random instance for the whole project. This Random can be initialized
 * in the Tools constructor with a fixed seed to allow reproducible results.
 * @author Roeland Voorrips
 */
public class Tools {

    public static double MAXMORGAN = 25.0;

    public Random rand;

    public Tools(long randomSeed) {
        rand = new Random(randomSeed);
        GAMMAFACTOR = GAMMAFACTOR_DEFAULT;
    }

    public int[] shuffle0123() {
        return shuffleArray(new int[] {0,1,2,3});
    }

    public int[] shuffleArray(int[] source) {
        int[] dest = new int[source.length];
        System.arraycopy(source, 0, dest, 0, source.length);
        //loop invariant: the tail of the sample array is randomized.
        //Initally the tail is empty; at each step move a random
        //element from front of array into the tail, then decrement boundary of tail
        int last = dest.length-1;   //last is maximal index of elements not in the tail

        while (last > 0){
            // Select random index in range 0..last, and swap its contents with those at last
            // The null swap is allowed; it should be possible that sample[k] does not change
            int ix = rand.nextInt(last+1);
            int tmp = dest[last];
            dest[last] = dest[ix];
            dest[ix] = tmp;
            last--;
        }
       return dest;
    } //shuffleArray

    /**
     * seq(n): 
     * @param n length of the result
     * @return an array of n integers from 0 to n-1
     */
    public static int[] seq(int n) {
        int[] result = new int[n];
        for (int i=0; i<n; i++) {
            result[i] = i;
        }
        return result;
    }
    
    /**
     * seq(from,to):
     * @param from integer, start of sequence (inclusive)
     * @param to integer, end of sequence (inclusive)
     * @return array with all integers from from to to, of length at least 1,
     * increasing or decreasing depending on to being larger or smaller than 
     * from
     */
    public static int[] seq(int from, int to) {
        int n = Math.abs(to-from)+1;
        int[] result = new int[n];
        int step = to<from ? -1 : 1;
        int m = from - step;
        for (int i=0; i<n; i++) {
            result[i] = m += step;
        }
        return result;
    }

    public static String printIntArray(int[] array) {
        if (array==null || array.length==0) {
            return "empty";
        } else {
            String s = "" + array[0];
            for (int i=1; i<array.length; i++) {
                s += "\t" + array[i];
            }
            return s;
        }
    }

    public static String printBooleanArray(boolean[] array) {
        if (array==null || array.length==0) {
            return "empty";
        } else {
            String s = "" + (array[0]?"1":"0");
            for (int i=1; i<array.length; i++) {
                s += array[i]?"\t1":"\t0";
            }
            return s;
        }
    }

    /**
     * haldaneToRecomb calculates the recombination fraction corresponding to
     * a given map distance according to the Haldane mapping function
     * Sybenga(1972) p.283
     * @param morgan the map distance
     * @return the recombination fraction
     */
    public static double haldaneToRecomb(double morgan) {
        return 0.5*(1.0-Math.exp(-2*Math.abs(morgan)));
    }

    /**
     * recombToHaldane is the inverse of haldaneToRecomb
     * @param recomb the recombination fraction
     * @return the map distance
     */
    public static double recombToHaldane(double recomb) {
        if (recomb<=0.0) return 0.0;
        if (recomb>=0.5) return Double.NaN;
        double x = 1.0 - 2*recomb;
        x = -0.5*Math.log(x);
        return x<MAXMORGAN ? x : MAXMORGAN;
    }

    /**
     * kosambiToRecomb calculates the recombination fraction corresponding to
     * a given map distance according to the Kosambi mapping function
     * Sybenga(1972) p.283
     * @param morgan the map distance
     * @return the recombination fraction
     */
    public static double kosambiToRecomb(double morgan) {
        return 0.5*Math.tanh(2*Math.abs(morgan));
    }

    /**
     * recombToKosambi is the inverse of kosambiToRecomb
     * @param recomb the recombination fraction
     * @return the map distance
     */
    public static double recombToKosambi(double recomb) {
        if (recomb<=0.0) return 0.0;
        if (recomb>=0.5) return Double.NaN;
        double x = 1.0 - 2*recomb;
        x = 0.25*Math.log((1.0+2*recomb)/x);
        return x<MAXMORGAN ? x : MAXMORGAN;
    }

    public static double mapdistToRecomb(double morgan, boolean chiasmaInterference) {
        if (chiasmaInterference) return kosambiToRecomb(morgan);
        else return haldaneToRecomb(morgan);
    }

    public static double recombToMapdist(double recomb, boolean chiasmaInterference) {
        if (chiasmaInterference) return recombToKosambi(recomb);
        else return recombToHaldane(recomb);
    }  
    
    /**
     * recombBivalentToRecombParallelQuadrivalent calculates the recombination
     * expected in a parallel quadrivalent at the distance where a bivalent
     * has the given input recombination.
     * According to Sved, Heredity 19: 585–596, 1964, formula(2) p. 592
     * @param bivalRec recombination fraction in a bivalent
     * @return corresponding recombination fraction in a parallel quadrivalent
     */
    public static double recombBivalentToRecombParallelQuadrivalent(double bivalRec) {
        return 0.75*(1-Math.pow((1-2*bivalRec),2.0/3.0));
    }
    
    /**
     * svedparallelToRecomb calculates the recombination fraction corresponding to
     * a given map distance according to the formula given by Sved (Heredity 19: 
     * 585–596, 1964) just above formula (2) p. 592
     * @param morgan the map distance
     * @return the recombination fraction
     */
    public static double svedparallelToRecomb(double morgan) {
        double m = morgan/4; //in a parallel quadrivalent there are 4 chiasmata per morgan
        return 0.75*(1-Math.exp(-2*m/3));
    }
    
    /**
     * recombToSvedparallel is the inverse of svedparallelToRecomb
     * @param recomb the recombination fraction
     * @return the map distance
     */
    public static double recombToSvedparallel(double recomb) {
        double m = -1.5*(Math.log(1-4*recomb/3));
        return m/4; //in a parallel quadrivalent there are 4 chiasmata per morgan
    }

    /**
     * recombBivalentToRecombCrossQuadrivalent calculates the recombination
     * expected in a cross-type quadrivalent at the distance where a bivalent
     * has the given input recombination.
     * According to Sved, Heredity 19: 585–596, 1964, formula(3) p. 594
     * Note that Sved assumes uniform distribution of the chromosome exchange
     * point which is different from the effect in our approach: we have a
     * unimodal distribution with the maximum at the chromosome center.
     * Sved's result is therefore not correct in our situation.
     * @param bivalRec recombination fraction in a bivalent
     * @return corresponding recombination fraction in a parallel quadrivalent
     */
    public static double recombBivalentToRecombCrossQuadrivalent(double bivalRec) {
        return bivalRec/(Math.log(1-2*bivalRec)) + bivalRec/2 + 0.5;
    }

    /**
     * svedcrosstypeToRecomb calculates the recombination fraction corresponding to
     * a given map distance according to the formula given by Sved (Heredity 19: 
     * 585–596, 1964) just before formula (3), end of p. 592
     * Note the remark at recombBivalentToRecombCrossQuadrivalent
     * @param morgan the map distance
     * @return the recombination fraction
     */
    public static double svedcrosstypeToRecomb(double morgan) {
        double m = morgan/2; //in a cross-type quadrivalent there are 4 chiasmata per morgan in each arm
        return 0.75 - 0.25*Math.exp(-m) - 0.5/m*(1-Math.exp(-m));
    }
    
    /**
     * ranExp gets a random number from the standard exponential distribution
     * @return
     */
    public double ranExp(double mean) {
        return  -mean*Math.log(1.0-rand.nextDouble());
    }

    /**
     * ranGamma returns a random sample from a Gamma distribution
     * (source: http://vyshemirsky.blogspot.com/2007/11/sample-from-gamma-distribution-in-java.html)
     * mean = k*theta, variance = k*theta*theta
     * Special cases:
     * k=1: exponential distribution with parameter 1/theta
     * theta=2: chi-square distribution with d.f. = 2*k
     * For getting the distance from one chiasma to the next with mean distance m,
     * taking interference into account:
     * k about 5 to 6 (McPeek&Speed 1995, human data) or
     * k about 2.63 (Kosambi mapping function; own simulations),
     * theta = m/k
     * For the first chiasma the exponential distribution with mean m should be used
     * (but in cross-type quadrivalents we start at both ends of each
     * chromosome; what to do here?) -->
     * Reasoning:
     * a) if no interference, any location beyond opposite chiasma is acceptable
     * b) the cumulative distribution with interference (with GAMMAFACTOR=2.63)
     *    is smaller than the one without interference (the exponential distribution)
     *    for distances under 0.1375 Morgan
     * c) so only in that region interference (rejection of new chiasma) should occur
     * d) we let the probability of rejection decrease quadratically from 1 at
     *    distance 0 to 0 at approx 0.1375, with the top of the parabola at
     *    that 0 point (approx x=0.1375)
     * e) all this is an approximation so we work with round parameters:
     *    P(rejection) = 64*x*x - 16*x + 1; 1 at x=0, 0 at x=0.125 Morgan
     *
     * @param k     : >0, shape parameter
     * @param theta : >0, scale parameter
     * @return
     */
    public double ranGamma(double k, double theta) {
        boolean accept = false;
        if (k < 1) {
             // Weibull algorithm
             double c = (1 / k);
             double d = ((1 - k) * Math.pow(k, (k / (1 - k))));
             double u, v, z, e, x;
             do {
              u = rand.nextDouble();
              v = rand.nextDouble();
              z = -Math.log(u);
              e = -Math.log(v);
              x = Math.pow(z, c);
              if ((z + e) >= (d + x)) {
               accept = true;
              }
             } while (!accept);
             return (x * theta);
        }
        else {
             // Cheng's algorithm
             double b = (k - Math.log(4));
             double c = (k + Math.sqrt(2 * k - 1));
             double lam = Math.sqrt(2 * k - 1);
             double cheng = (1 + Math.log(4.5));
             double u, v, x, y, z, r;
             do {
              u = rand.nextDouble();
              v = rand.nextDouble();
              y = ((1 / lam) * Math.log(v / (1 - v)));
              x = (k * Math.exp(y));
              z = (u * v * v);
              r = (b + (c * y) - x);
              if ((r >= ((4.5 * z) - cheng)) ||
                                (r >= Math.log(z))) {
               accept = true;
              }
             } while (!accept);
             return (x * theta);
        }
    }

    static double GAMMAFACTOR_DEFAULT = 2.63; //2.63 gives a good fit of Kosambi mapping function
    /* with GAMMAFACTOR = 1.0 the non-interference/Haldane situation is obtained
     * but this in implemented using the (equivalent but faster) exponential distribution
     */
    double GAMMAFACTOR;
   
    /**
     * ranDistInterference
     * @param meandist the desired average distance between chiasmata
     * @return a random sample from a gamma distribution, to be used as
     * the distance from one chiasma to the next. This will lead to a
     * good correspondence with the Kosambi mapping function.
     */public double ranDistInterference(double meandist) {
        return ranGamma(GAMMAFACTOR, meandist/GAMMAFACTOR);
    }

    private static double mdiff (double L, double m) {
        //used by calcChiasmaDist
        return Math.abs(m - 0.5 / (1-Math.exp(-L/m)));
    }

    /**
     * calcChiasmaDist calculates the necessary average distance of chiasmata
     * m, such that there are on average 2*L chiasmata on a chromosome with
     * length L, IF all bivalents without chiasmata are rejected.
     * In that case, m == 0.5 / (1-exp(-L/m))
     * NOTE that this will not lead to a correspondence with the Haldane or
     * Kosambi mapping functions; in particular the total map length of
     * the chromosome in a bivalent will always be infinite
     * @param chromlength
     * @return
     */
    public static double calcChiasmaDist(double chromLength) {
         if (chromLength<=0.500) {
             return Double.NaN;
         }
         if (chromLength>5.0) {
             return 0.5;
         }
         double toler = 0.00001;
         double lastdiff = mdiff(chromLength,0.5);
         //first find initial minimum:
         double m1 = 10.5;
         while (mdiff(chromLength,m1)<lastdiff) {
             lastdiff = mdiff(chromLength,m1);
             m1 += 10.0;
         }
         double m2,m3;
         if (m1==10.5) {
             m2 = 5.5;
             m3 = 0.5;
         } else {
             m2 = m1 - 10.0;
             m3 = m1 - 20.0;
        }
         while (m1-m3>toler) {
             if (mdiff(chromLength,m1-0.5*(m1-m2))<mdiff(chromLength,m2)) {
                 m3 = m2;
                 m2 = m1-0.5*(m1-m2);
             }
             else if (mdiff(chromLength,m3+0.5*(m2-m3))<mdiff(chromLength,m2)) {
                 m1 = m2;
                 m2 = m3+0.5*(m2-m3);
             } else {
                 m1 = m1-0.5*(m1-m2);
                 m3 = m3+0.5*(m2-m3);
             }
         }
         return m2;
    } //calcChiasmaDist
    
    public static double sampleStDev(double sum, double sumSquares, int nobs) {
        return Math.sqrt((sumSquares-sum*sum/nobs)/(nobs-1));
    } //sampleStDev}

    public static String[] readWords(String line) {
        if (line==null || line.trim().isEmpty()) {
            return new String[0]; //not null: allows result to be tested for length
        }
        return line.trim().split("\\s+"); //split on all whitespace sequences
        /* Problem with split on an empty String:
         * returns a String[1] with element[0]==""
         * Complicates testing for emptiness
         */
    } //readWords

    public static ArrayList<String> readWordsQ(String line, char quote) {
        //based on readWoord, reads quoted parts as one word including spaces
        ArrayList<String> words = new ArrayList<String>();
        int p=0, q, le=line.length();
        while (p<le) {
            while (p<le && line.charAt(p)<=' ') p++;
            if (p<le && line.charAt(p)==quote) {
                //read quoted word until next quote or EOLN
                // (may contain any characters including blanks, except quote)
                q=++p; //q and p past start-quote
                while (p<le && line.charAt(p)!=quote) p++;
                if (p>q) {
                    words.add(line.substring(q, p));
                }
                p++; //past end-quote
            }
            else { //not quoted
                q = p++;
                while (p<le && line.charAt(p)>' ') p++;
                if (p>q) {
                    words.add(line.substring(q, p));
                }
            }
        }
        return words;
    } //readWordsQ

}
