/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package PedigreeSim;

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
    //(the exchange interval is the interval between the last chiasma in the start arms
    //and the last chiasma in the end arms of the quadrivalent):
    static int freqTableLength = 20;
    int crossQuadrivalentCount;
    int paralQuadrivalentCount;
    int bivalentCount;
    int noExchangeLimCount; /*counts arm-based quadrivalents where neither the 
     *      start nor the end if the exchange interval is defined (no chiasmata)*/
    int oneExchangeLimCount;/* counts the arm-based quadrivalents where only the
     *      start or the end of the exchange interval is defined (chiasmata only
     *      in the start arms or only in the end arms*/
    int twoExchangeLimCount;/* counts the arm-based quadrivalents where both the
     *      start and the end of the exchange interval are defined (chiasmata both
     *      in the start arms and the end arms*/
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
        this.prefPairingProb = prefPairingProb;
        this.fracQuadrivalents = fracQuadrivalents;
    }

    public double getFracQuadrivalents() {
        return fracQuadrivalents;
    }

    public double getPrefPairingProb() {
        return prefPairingProb;
    }


}
