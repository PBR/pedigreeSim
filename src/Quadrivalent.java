/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package PedigreeSim;

/**
 * A Quadrivalent is composed of four homologous chromosomes (8 chromatids).
 * @author Roeland Voorrips
 */
public class Quadrivalent extends Multivalent {

    private static int[][] arm = new int[][] {
        {0,1,6,7}, //chromatids with start in arm 0
        {0,1,2,3}, //chromatids with end in arm 1
        {2,3,4,5}, //chromatids with start in arm 2
        {4,5,6,7} //chromatids with end in arm 3
    }; //see discussion in Javadoc of doMeiosis

    /**
     * The order of the chromosomes is NOT randomized in the constructor!
     * Instead it is determined by the caller (Individual.doMeiosis)
     * because it depends on the preferential pairing of chromosome
     * ends, which is known in the Individual.
     * @param haplostruct
     * @throws Exception
     */
    public Quadrivalent(HaploStruct[] haplostruct) throws Exception {
        super();
        boolean ok = haplostruct!=null && haplostruct.length==4 &&
                haplostruct[0]!=null && haplostruct[0].getChrom()!=null;
        for (int h=1; h<4; h++) {
            ok = ok && haplostruct[h]!=null &&
                    haplostruct[0].getChrom().equals(haplostruct[h].getChrom());
        }
        if (!ok) {
            throw new Exception("Error in Quadrivalent constructor parameters");
        }
        this.chrom = haplostruct[0].getChrom();
        this.popdata = chrom.getPopdata();
        this.tools = popdata.tools;
        this.rand = tools.rand;
        this.haplostruct = haplostruct;
    }

    /**
     * getStartArm:
     * gives the arm number in which the start of the given chromatid lies
     * @param chromatid
     * @return 
     */
    private int getStartArm(int chromatid) {
        return (chromatid>=2 && chromatid<=5) ? 2 : 0;
    }

    /**
     * getEndArm:
     * gives the arm number in which the end of the given chromatid lies
     * @param chromatid
     * @return 
     */
    private int getEndArm(int chromatid) {
        return (chromatid<=3) ? 1 : 3;
    }

    //some global variables for use in doCrossingOver, generateChiasma, getCentromereSortOrder:
    private HaploStruct[] slots;
    private int side; //0=start: arm 0 and 2, 1=end: arm 1 and 3
    private double sideFact; //multiplication factor to reverse directions
    private double mindist; //dist to nearest of two opposing chiasmata
    private int chiasmaCount;
    private double pos;       //position of the new chiasma
    private int[] chromatid = new int[2];   //the two slots involved in the new chiasma
    private double[] exchangeLim = new double[2]; //start (0) and end (1) limits of the exchange interval
    private double[] lastpos = new double[4]; //position of last chiasma in each arm
    private boolean[] armFinished = new boolean[4];
    private boolean crossoverFinished; //true if whole process finished

    /**
     * generateChiasma: internal function
     * Generates a new chiasma in the specified arm, checks for conflicts
     * and adjusts the remaining exchange interval according to
     * doCrossingOver with quadriMethod=2
     * @param a the number of the arm (0..3) where the chiasma will be generated
     * @param chiasmaDist the average distance between two chiasmata (either
     * CHIASMADIST if popdata.allowNoChiasmata==true, or else
     * Tools.calcCiasmaDist)
     * The result affects the private global variables above:
     * this.pos is the position of the new chiasma
     * this.chromatid[0..1] indicate the two slots involved in the new chiasma
     * this.armFinished[a] tells if the new chiasma is valid
     * this.lastpos[a] is updated if the chiasma is valid
     * this.exchangeLim[side] is updated if needed
     */
    private void generateChiasma(int a, double chiasmaDist) {
        assert(!crossoverFinished);
        assert(!armFinished[a]);
        assert(slots.length==8);
        //select two of the four slots (chromatids) in this arm for chiasma formation:
        chromatid[0] = arm[a][rand.nextInt(2)];
        chromatid[1] = arm[a][2+rand.nextInt(2)];
        side = a % 2; //0=start in arm 0 and 2, 1=end in arm 1 and 3
        int oppoArm1 = 1-side; //the two opposing arms: 1 and 3, or 0 and 2
        int oppoArm2 = 3-side;
        /* note that:
         * (1-side) is the opposing side (end=1 vs start=0)
         * the other arm of the same side is (a+2)%4
         */
        sideFact = side==0 ? 1.0 : -1.0; //multiplication factor to reverse directions
        //calculate position of new chiasma:
        if (Double.isNaN(lastpos[a])) {
            if (side==0) pos = chrom.getStartPos() + distToFirstChiasma(chiasmaDist);
            else pos = chrom.getEndPos() - distToFirstChiasma(chiasmaDist);
        }
        else {
            pos = lastpos[a] + sideFact*distToNextChiasma(chiasmaDist);
        }
        mindist = Double.NaN; //dist to nearest of two opposing chiasmata
        if (!Double.isNaN(lastpos[oppoArm1])) {
            mindist = sideFact*(lastpos[oppoArm1]-pos);
        }    
        if (!Double.isNaN(lastpos[oppoArm2]) &&
            (Double.isNaN(mindist) || mindist>sideFact*(lastpos[oppoArm2]-pos))) {
                mindist = sideFact*(lastpos[oppoArm2]-pos);
        }
    } //generateChiasma
    
    private void testChiasma2(int a) {
        //NOTE: must be called after generateChiasma; 
        //there side, sideFact, pos, mindist are calculated
        double tempLim = exchangeLim[side]; //temporary exchangeLim at this side
        //first: check if rejected because of opposing chiasma
        //(possibly with interference):
        if (!Double.isNaN(mindist)) {
            //there are opposing chiasmata; conflict?
            if (rejectDist(mindist)) {
                // conflict with opposing chiasma
                armFinished[a] = true;
                if (mindist<=0) {
                    //situation 2b: new chiasma generated beyond opposing chiasma
                    tempLim = pos + sideFact*mindist; //at opposing chiasma;
                }
                else {
                    //situation 2c: new chiasma generated at this side of 
                    //opposing chiasma but rejected because of interference
                    tempLim = pos;
                }
            }
            else {
                //opposing chiasmata but no conflict:
                tempLim = pos;
            }
        }
        
        //check if pos not beyond opposite exchangeLim (situation 2a or 2d):
        if (sideFact*(pos-exchangeLim[1-side])>=0) {
            armFinished[a] = true;
            tempLim = exchangeLim[1-side];
        }
        
        if (armFinished[a]) {
            crossoverFinished = true;
            for (int i=0; i<4; i++) {
                crossoverFinished = crossoverFinished && armFinished[i];
            }
            //false if at least one armFinished false
        }

        //update exchangeLim on this side if tempLim beyond it:
        if (sideFact*(tempLim-exchangeLim[side])>0) {
            exchangeLim[side] = tempLim;
        }
        
    } //testChiasma2

    private void testChiasma13(int a) {
        //NOTE: must be called after generateChiasma; 
        //there side, sideFact, pos, mindist are calculated
        //first: check if rejected because of opposing chiasma
        //(possibly with interference):
        if (!Double.isNaN(mindist) && rejectDist(mindist)) {
            // conflict with opposing chiasma
            armFinished[a] = true;
        }
        else { //test if chiasma beyond opposite chromosome end:
            armFinished[a] = (side==0 && pos>=chrom.getEndPos()) ||
                             (side==1 && pos<=chrom.getStartPos());
        }    
        
        if (armFinished[a]) {
            crossoverFinished = true;
            if (popdata.quadriMethod==1) {
                for (int i=0; i<4; i++) {
                    crossoverFinished = crossoverFinished && armFinished[i];
                }
                //false if at least one armFinished false
            } 
            //else quadriMethod==3: crossoverFinished as soon as one arm finished
        }
        else { //!armFinished[a]
            //update exchangeLim on this side if pos beyond it:
            if (sideFact*(pos-exchangeLim[side])>0) {
                exchangeLim[side] = pos;
            }
        }
    } //testChiasma13

    private void makeChiasma(int a) {
        assert(!armFinished[a]);
        chiasmaCount++;
        //insert the recombination into both slots:
        slots[chromatid[0]].insertSegment(pos, chromatid[1]);
        slots[chromatid[1]].insertSegment(pos, chromatid[0]);
        //update lastpos:
        lastpos[a] = pos;
    } //makeChiasma 
        

    /**
     * doCrossingOver performs crossing-over for a cross-type Quadrivalent.
     * First we define 8 new HaploStructs ("slots") to represent the chromatids,
     * each with just one segment of one of four pseudo-founderalleles.
     * The four chromosomes are each composed of 2 chromatids ("slots"),
     * 0+1, 2+3, 4+5 and 6+7.
     * Picture the Quadrivalent as a cross, with a left, up, right and down arm.
     * In each arm two of the chromosome ends are paired as follows:
     * left arm: start of slots 0,1,6,7
     * up arm: end of slots 0,1,2,3
     * right arm: start of slots 2,3,4,5
     * down arm: end of slots 4,5,6,7
     * (in the following we say that the up and down arms are the "opposing
     * arms" of the left and right arms, and v.v.)
     * We define 3 doCrossingOver methods for a cross-type Quadrivalent:
     * 1. Chiasmata are generated from the ends of each arm, until they
     * conflict with an existing chiasma in an opposite arm or are located
     * beyond the opposite end. When a conflicting chiasma has occurred
     * no new chiasmata are generated for that arm, but in the other arms
     * the process continues, until in all arms a conflicting chiasma has 
     * occurred.
     * 2. Chiasmata are generated from the ends of each arm until the first
     * conflict occurs. A conflict can occur because:
     * 2a. the new chiasma is beyond the opposite end of the chromosome
     * (this is only possible if no chiasmata have occurred in the two opposing 
     * arms, else situation 2b applies).
     * In this case in effect there are two bivalents, and for one of these
     * bivalents the crossing-over process has just finished.
     * Implementation: if the arm has the chromosome start, then the
     * exhangeStart becomes chromosomeEnd, exchangeEnd was already chromosomeEnd
     * and armFinished is true for this arm. Valid chiasmata can only be 
     * generated in the other arm that has also the chromosome start; until 
     * here also a chiasma is generated beyond the chromosomeEnd.
     * 2b. the new chiasma is beyond the nearest chiasma in one of the two
     * opposing arms. In that case the chromosome exchange point is now fixed
     * at that chiasma position and no new chiasmata can be generated in the
     * two arms involved in the conflict. In the other two arms further chiasmata
     * are generated until they conflict with the exchange point.
     * Implementation: if the current arm has the chromosome start, the
     * exchangeStart and exchangeEnd are now set to the position of the 
     * conflicting chiasma and no more chiasmata can be generated in these
     * two arms.
     * 2c. the new chiasma is generated on before the nearest chiasma on the
     * two opposing arm, but close to it and due to interference it is conflicting
     * (this can only occur when popdata.chiasmaInterference is true). In that
     * case no new chiasmata are generated for this arm, but the failed
     * chiasma defines the new start point of the chromosome exchange interval.
     * Implementation: if the current arm has the chromosome start,
     * exchangeStart is set to the failed position and armFinished is
     * true.
     * 2d. the new chiasma is generated beyond the position of a previous
     * conflicting chiasma as discussed in 2c, but not conflicting with
     * an actual chiasma. Then the start of the exchange interval is
     * moved to that position, which is also the end of the interval
     * Implementation: if the current arm has the chromosome start, then 
     * this situation occurs if no conflicts with opposing chiasmata occur
     * but the new chiasma is beyond exchangeEnd. Now exchangeStart is set
     * to exchangeEnd and armFinished is true. This is the same as what
     * happens in situation 2a.
     * The reasoning behind method 2 (a, b and c) is that in order to generate
     * a chiasma, even if it conflicts, the pairing of the chromosomes in that
     * arm must already have been completed up to that point or up to the 
     * chromosome exchange point.
     * QUESTION: does the chromosome exchange point cause interference like
     * a chiasma? 
     * ANSWER: that depends on whether the interference is mainly caused by the
     * mechanical bending stress (then: yes) or mainly by the break/repair
     * process (then: no). But a practical consideration is that it is often
     * not precisely known where the exchange point is, so interference
     * cannot be modelled.
     * 3. Chiasmata are generated from the ends of each arm until the first
     * conflict occurs. When this happens the crossing-over process is
     * ended for the complete quadrivalent.
     * This seems to be too stringent: the amount of recombination at the
     * ends of the arms is about half that in a bivalent.     
     * @return a HaploStruct[] with 8 elements: the 8 recombined chromatids,
     * (the slots) still in the original order.
     * @throws Exception
     */
    @Override
    protected HaploStruct[] doCrossingOver() throws Exception {
        boolean paralMv = popdata.parallelMultivalents==0.0 ? false :
                          popdata.parallelMultivalents==1.0 ? true :
                rand.nextDouble()<popdata.parallelMultivalents;
        if (paralMv) {
            return super.doCrossingOver(); //parallel quadrivalent
        }
        
        //cross-type quadrivalent :
        do {
            double chiasmaDist = CHIASMADIST;
            if (!popdata.allowNoChiasmata) {
                chiasmaDist = Tools.calcChiasmaDist(chrom.getLength());
            }
            
            slots = new HaploStruct[8];
            for (int s=0; s<8; s++) {
                slots[s] = new HaploStruct(chrom, -1); //the first segment will be ignored
            }
            for (int a=0; a<4; a++) {
                lastpos[a] = Double.NaN;
                //if (a==0 || a==2) exchangelim[a] = chrom.getStartPos();
                //else exchangelim[a] = chrom.getEndPos();
                armFinished[a] = false;
            }
            crossoverFinished = false;
            chiasmaCount = 0;
            exchangeLim[0] = chrom.getStartPos();
            exchangeLim[1] = chrom.getEndPos();
            if (popdata.quadriEachArm) {
                int[] aix; //the sequence of the arms in each cycle
                do {
                    aix = tools.shuffle0123();
                    for (int a=0; a<4; a++) {
                        if (!armFinished[aix[a]]) {
                            generateChiasma(aix[a],chiasmaDist);
                            if (popdata.quadriMethod==2) testChiasma2(aix[a]);
                            else testChiasma13(aix[a]);
                            if (!armFinished[aix[a]]) makeChiasma(aix[a]);
                        } // !armFinished
                    } //for a
                } while (!crossoverFinished); 
            }
            else { //not popdata.quadriEachArm
                //(same as above without aix and the for a - loop)
                do {
                    int a = rand.nextInt(4); //which arm for the next chiasma
                        if (!armFinished[a]) {
                            generateChiasma(a, chiasmaDist);
                            if (popdata.quadriMethod==2) testChiasma2(a);
                            else testChiasma13(a);
                            if (!armFinished[a]) makeChiasma(a);
                        } // !armFinished
                } while (!crossoverFinished);
            }
        } while (!(popdata.allowNoChiasmata || allChromChiasmata(slots)));
        
        if (popdata.testMode) {
            ((TetraploidChromosome) chrom).crossQuadrivalentCount++;
            ((TetraploidChromosome) chrom).quadChiasmaSum += chiasmaCount;
            ((TetraploidChromosome) chrom).quadChiasmaSS += chiasmaCount*chiasmaCount;
            //statistics on midpoint and length of exchange interval:
            if (chiasmaCount==0 ||
                (exchangeLim[0]==chrom.getStartPos() && 
                 exchangeLim[1]==chrom.getEndPos())) {
                ((TetraploidChromosome) chrom).noExchangeLimCount++;
            } else if (exchangeLim[0]==chrom.getStartPos() || 
                       exchangeLim[1]==chrom.getEndPos()) {
                ((TetraploidChromosome) chrom).oneExchangeLimCount++;
            }
            else {
                ((TetraploidChromosome) chrom).twoExchangeLimCount++;
                double x = (exchangeLim[0]+exchangeLim[1]) / 2.0; //midpoint of exchange interval
                ((TetraploidChromosome) chrom).exchangeMidSum += x;
                ((TetraploidChromosome) chrom).exchangeMidSS += x*x;
                int i = (int)(Math.floor((x-chrom.getStartPos())*TetraploidChromosome.freqTableLength 
                        / chrom.getLength()));
                ((TetraploidChromosome) chrom).exchangeMidFreq[i]++;
                x = exchangeLim[1]-exchangeLim[0]; //length of exchange interval
                ((TetraploidChromosome) chrom).exchangeLengthSum += x;
                ((TetraploidChromosome) chrom).exchangeLengthSS += x*x;
                i = (int)(Math.floor((x-chrom.getStartPos())*TetraploidChromosome.freqTableLength 
                        / chrom.getLength()));
                ((TetraploidChromosome) chrom).exchangeLengthFreq[i]++;
            }    
        }
        return slots;
    } //doCrossingOver
        
    /**
     * Quadrivalent.getCentromereSortOrder:
     * see for a full explanation Multivalent.getCentromereSortOrder
     * @return 
     */
    @Override
    int[] getCentromereSortOrder() {
        boolean paired = popdata.pairedCentromeres==0.0 ? false :
                popdata.pairedCentromeres==1.0 ? true :
                rand.nextDouble()<popdata.pairedCentromeres;
        if (!paired) {
            /* this assumes that all centromere pairs have equal chance
             * of being formed at first meiotic division
             */
            return super.getCentromereSortOrder();
        }
        /* alternatively, this assumes that of each set of paired centromeres
         * one goes to each opposite pole at the first meiotic division
         */
        int[][] gam; //which centromeres in which gamete
        if (rand.nextBoolean()) {
            //Note that exchangeLim has been actualized in testChiasma
            boolean centromeresInStartArm = chrom.getCentromerePos() <= exchangeLim[1]; 
            //adjacent configuration:
            if (centromeresInStartArm) {
                if (rand.nextBoolean()) 
                     gam = new int[][] { {0,1}, {0,1}, {2,3}, {2,3} };
                else gam = new int[][] { {2,3}, {2,3}, {0,1}, {0,1} };
            }
            else { //centromeres in end arm
                if (rand.nextBoolean())
                     gam = new int[][] { {0,3}, {0,3}, {1,2}, {1,2} };
                else gam = new int[][] { {1,2}, {1,2}, {0,3}, {0,3} };
            }
        }
        else {
            //alternate configuration, independent of arm where centromeres are:
            if (rand.nextBoolean())
                 gam = new int[][] { {0,2}, {0,2}, {1,3}, {1,3} };
            else gam = new int[][] { {1,3}, {1,3}, {0,2}, {0,2} };
        }
        // randomize the order in each gamete:
        for (int g=0; g<4; g++) {
            gam[g] = tools.shuffleArray(gam[g]);
        }
            //next, put gam[0..3] consecutively into centro:
        int[] centro = new int[2*haplostruct.length];
        for (int g=0; g<4; g++) {
            System.arraycopy(gam[g], 0, centro, g*(centro.length/4), centro.length/4);
        }
        return centro;
    } //getCentromereSortOrder

}
