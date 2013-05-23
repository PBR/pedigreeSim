/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package PedigreeSim;

import JSci.maths.statistics.BinomialDistribution;

/**
 * TetraploidChromosome extends Chromosome with two fields:
 * prefPairingProb and fracQuadrivalents
 * @author Roeland Voorrips
 */
public class TetraploidChromosome extends Chromosome {
    
    /**
     * prefPairingProb: 
     * probability of preferential pairing vs random pairing of the
     * telomeres
     */
    private double prefPairingProb;
    
    /**
     * fracQuadrivalents:
     * fraction of meioses in which the four homologues form one
     * Quadrivalent instead of two Bivalents. 
     * Only used if popdata.naturalPairing==false (but it must be specified
     * even if it is not used)
     */
    private double fracQuadrivalents; 

    //Arrays for use in test mode:
    //arrays for storing the cumulative counts of homozygosity for founder alleles at each marker, 
    //index=markers
    int[] founderHomAll; //over all gametes
    int[] founderHomSel; //for one random gamete per meiosis
    //int[] founderHomCum; //founderHomSel cumulated over all iterations
    double[] founderHomFrCum; //freq of homozygozity cumulated over all ietration
    double[] founderHomFrSqCum; //squared frq of homoz cumulated over all iterations
    //variables/arraysfor storing a frequency distribution and sums/sums-of-squares
    //of the chrosomome exchange interval in quadrivalents
    //(the exchange interval is the interval between the last chiasma in the head arms
    //and the last chiasma in the tail arms of the quadrivalent):
    static int freqTableLength = 20;
    int[] quadrivalentConfigCount; //counts meioses with 0, 1, ... quadrivalents
    int crossQuadrivalentCount; //counts all cross quadrivalents
    int paralQuadrivalentCount; //counts all parallele quadrivalents
    int bivalentCount; //counts all bivalents
    int noExchangeLimCount; /*counts arm-based quadrivalents where neither the 
     *      start nor the end if the exchange interval is defined (no chiasmata)*/
    int oneExchangeLimCount;/* counts the arm-based quadrivalents where only the
     *      start or the end of the exchange interval is defined (chiasmata only
     *      in the head arms or only in the tail arms*/
    int twoExchangeLimCount;/* counts the arm-based quadrivalents where both the
     *      start and the end of the exchange interval are defined (chiasmata both
     *      in the head arms and the tail arms*/
    int[] exchangeMidFreq; //freq.table of the position of the centre of the exchange interval
    int[] exchangeLengthFreq; //freq.table of the length of the exchange interval
    double exchangeMidSum; //sum of the centre positions of the exchange intervals
    double exchangeMidSS; //sum of squares
    double exchangeLengthSum; //sum of the lengths of the exchange intervals
    double exchangeLengthSS; //sum of squares
    //variables for storing the sums/sums-of-squares for the number of chiasmata
    //in quadrivalents:
    int quadChiasmaSum; //sum
    int quadChiasmaSS; //sum of squares

    public TetraploidChromosome(String chromName, double length, double centromerePos,
            double prefPairingProb, double fracQuadrivalents,
            PopulationData popdata) {
        super(chromName, length, centromerePos, popdata);
        assert(popdata.ploidy > 2);
        this.prefPairingProb = prefPairingProb;
        this.fracQuadrivalents = fracQuadrivalents;
        this.quadrivalentConfigCount = new int[1 + popdata.ploidy/4];
        setCumulBinomial();
    }

    public double getFracQuadrivalents() {
        return fracQuadrivalents;
    }

    public double getPrefPairingProb() {
        return prefPairingProb;
    }

    // cumulBinomial and setCumulBinomial are needed for ranQuadri:
    private double[] cumulBinomial;

    private void setCumulBinomial() {
        BinomialDistribution bindist =
                new BinomialDistribution(popdata.ploidy/4, fracQuadrivalents);
        cumulBinomial = new double[popdata.ploidy/4];
        for (int i=0; i<cumulBinomial.length; i++) {
            cumulBinomial[i]=bindist.cumulative(i);
        }
        //TODO remove debug code:
        //System.out.println("chrom "+this.getChromName()+"\tfracQuadri="+fracQuadrivalents);
        //for (int i=0; i<cumulBinomial.length; i++) System.out.println(i+"\t"+cumulBinomial[i]);
    }

    /**
     * ranQuadri determines a random number of Quadrivalents to be generated
     * given the ploidy of the population and fracQuadrivalents.
     * This kludge is needed because BinomialDistribution.inverse is flawed
     * (and it may be faster too for realistic ploidy values).
     * @return a random number (int) taken from a binomial distribution
     * with n = max. possible number of quadrivalents (=ploidy/4)
     * and p = fracQuadrivalents
     */
    int ranQuadri() {
        double p = tools.rand.nextDouble();
        int i = 0;
        while (i<cumulBinomial.length && cumulBinomial[i]<p) i++;
        return i;
    }

}
