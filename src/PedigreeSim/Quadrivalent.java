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

    private static final int[][] arm = new int[][] {
        {0,1,6,7}, //chromatids with head in arm 0
        {0,1,2,3}, //chromatids with tail in arm 1
        {2,3,4,5}, //chromatids with head in arm 2
        {4,5,6,7} //chromatids with tail in arm 3
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
     * getHeadArm:
     * gives the arm number in which the head of the given chromatid lies
     * @param chromatid
     * @return 
     */
    private int getHeadArm(int chromatid) {
        return (chromatid>=2 && chromatid<=5) ? 2 : 0;
    }

    /**
     * getTailArm:
     * gives the arm number in which the tail of the given chromatid lies
     * @param chromatid
     * @return 
     */
    private int getTailArm(int chromatid) {
        return (chromatid<=3) ? 1 : 3;
    }

    //some global variables for use in doCrossingOver, generateRecombination, getCentromereSortOrder:
    private HaploStruct[] slots;
    private int side; //0=head: arm 0 and 2, 1=tail: arm 1 and 3
    private int nearestArm; //opposite arm (0..3) with nearest recombination, -1=none
    private static final int[][] oppoArms = new int[][] {{1,3},{0,2}};
    private double sideFact; //multiplication factor to reverse directions
    private double mindist; //dist to nearest of two opposing recombinations
    private int recombCount;
    private double pos;       //position of the new recombination
    private final int[] chromatid = new int[2];   //the two slots involved in the new recombination
    private final double[] exchangeLim = new double[2]; //start (0) and end (1) limits of the exchange interval
    private final double[] lastpos = new double[4]; //position of last recombination in each arm
    private final boolean[] armFinished = new boolean[4];
    private boolean crossoverFinished; //true if whole process finished

    /**
     * generateRecombination: internal function
     * Generates a new recombination in the specified arm, checks for conflicts
     * and adjusts the remaining exchange interval according to
     * doCrossingOver with quadriMethod=2
     * @param a the number of the arm (0..3) where the recombination will be generated
     * @param recombDist the average distance between two recombinations (either
     * RECOMBDIST if popdata.allowNoRecomb==true, or else Tools.calcRecombDist)
     * The result affects the private global variables above:
     * this.pos is the position of the new recombination
     * this.side = 0 or 1 depending on the arm a
     * this.sidefact = 1 or -1 if side is 0 or 1
     * this.chromatid[0..1] indicate the two slots involved in the new recombination
     * this.mindist = distance of new recombination to nearest opposite recombination
     * or NaN if none
     * this.nearestArm is the opposite arm where the nearest recombination is, 
     * or -1 if none
     * The following private global variables are not changed in generateRecombination
     * but in either testRecomb2 or testRecomb13 which are called immediately afterwards:
     * this.armFinished[a] is set to true if the new recombination is generated 
     * beyond the chromosome end, beyond the nearest recombination or too close
     * to the nearest recombination in an opposite arm);
     * in testRecomb2 armFinished of the opposing arms are also affected
     * this.lastpos[a] is updated if the recombination is valid
     * this.exchangeLim[side] is updated if needed
     */
    private void generateRecombination(int a, double recombDist) {
        assert(!crossoverFinished);
        assert(!armFinished[a]);
        assert(slots.length==8);
        //select two of the four slots (chromatids) in this arm for recombination:
        chromatid[0] = arm[a][rand.nextInt(2)];
        chromatid[1] = arm[a][2+rand.nextInt(2)];
        side = a % 2; //0=head in arm 0 and 2, 1=tail in arm 1 and 3
        /* note that:
         * (1-side) is the opposing side (tail=1 vs head=0)
         * the other arm of the same side is (a+2)%4
         */
        sideFact = side==0 ? 1.0 : -1.0; //multiplication factor to reverse directions
        //calculate position of new recombination:
        if (Double.isNaN(lastpos[a])) {
            pos = chrom.getSidePos(side) + sideFact*distToFirstRecomb(recombDist);
        }
        else {
            pos = lastpos[a] + sideFact*distToNextRecomb(recombDist);
        }
        mindist = Double.NaN; //dist to nearest of two opposing recombinations
        nearestArm = -1; //none yet
        if (!Double.isNaN(lastpos[oppoArms[side][0]])) {
            mindist = sideFact*(lastpos[oppoArms[side][0]]-pos);
            nearestArm = oppoArms[side][0];
        }    
        if (!Double.isNaN(lastpos[oppoArms[side][1]]) &&
            (Double.isNaN(mindist) || mindist>sideFact*(lastpos[oppoArms[side][1]]-pos))) {
                mindist = sideFact*(lastpos[oppoArms[side][1]]-pos);
                nearestArm = oppoArms[side][1];
        }
    } //generateRecombination
    
    private void testRecomb2(int a) {
        //NOTE: must be called after generateRecombination; 
        //there side, sideFact, pos, mindist, nearestArm are calculated
        double tempLim = exchangeLim[side]; //temporary exchangeLim at this side
        double distExchangelim = sideFact * (pos-exchangeLim[1-side]); 
            //distance to opposite side of exchange interval: >0 if beyond, <0 if before 
        if (distExchangelim>=0.0 || 
                (!Double.isNaN(mindist) && rejectDist(mindist)) ) {
            /* new recombination pos conflicts with earlier recombinations or 
             * attempted recombinations or chromosome end on opposite branches
             * because one of following reasons:
             * 1. no opposite recombinations, new recombination beyond chrom end: 
             *    distExchangelim>=0
             * 2. attempted and/or actual opposite recombination(s),
             *    new recombination beyond closest position of these: 
             *    distExchangelim>=0
             * 3. actual opposite recombination(s), new recombination before these
             *    but rejected due to interference: distExchangelim<0
             *
             * Is 2 correct? If exchangeLim[1-side] was moved due to
             * a failed (interference) recombination, is a new recombination
             * beyond this limit then wrong?
             * Yes, because the failed recombination implied that the pairing
             * in the opposite arm had progressed to that point, so a new
             * recombination from the current arm at or beyond that position
             * is not possible any more
             */
                    
                            
            armFinished[a] = true;
            if (distExchangelim>=0.0) {
                /* new recombination beyond opposite side of exchange interval,
                 * reason 1 or 2 above:
                 * no new recombination, but exchange interval now closed
                 * (by moving this side to the opposite side)
                 */
                tempLim = exchangeLim[1-side];
                /* if tempLim now at opposite end of chrom the opposite arms 
                 * are finished, 
                 * or if not, the opposite arms with a recombination at 
                 * exchangeLim[1-side] are now also finished:
                 */
                for (int oa=0; oa<2;oa++) {
                    double lpoa;
                    if ( tempLim == chrom.getSidePos(1-side) ||
                         ( !Double.isNaN(lpoa=lastpos[oppoArms[side][oa]]) &&
                           (sideFact * (lpoa-tempLim)) <= 0 ) ) {
                        armFinished[oppoArms[side][oa]] = true;
                    }    
                }  
            } else {
                /* distExchangelim<0.0	
                 * only possible due to interference, reason 3 above:
                 * no new recombination, but reduce exchange interval by moving
                 * this side to pos of this attempted recombination (if that is
                 * beyond the current start of the interval)
                 */
                assert(popdata.chiasmaInterference);
                if (sideFact * (pos-tempLim) > 0) {
                    tempLim = pos; //truncatedPos;
                    //tempLim2Otherside = true;
                }
                /* the arm with the recombination which the interference took place 
                 * is now also finished:
                 */
                assert(nearestArm>-1);
                armFinished[nearestArm] = true;
            } 
            //this arm is now finished; check if all arms now finished:
            crossoverFinished = true;
            for (int i=0; i<4; i++) {
                crossoverFinished = crossoverFinished && armFinished[i];
            }
            //false if at least one armFinished false
        }
        else {
            /* distExchangelim<0.0 && 
             * (Double.isNaN(mindist) || !rejectDist(mindist)), meaning:
             * new recombination not beyond opposite end of exchange interval,
             * i.e. not beyond chromosome end and not beyond opposing recombination,
             * so on an allowed part of chromosome; AND
             * either no opposing recombination (Double.isNaN(mindist))
             * or no interference with existing opposing recombination (!rejectDist(mindist))
             * Therefore: no conflict, accept new recombination
             */
            lastpos[a] = pos;
            /* reduce exchange interval by moving
             * this side to pos of this recombination (if that is
             * beyond the current start of the interval):
             */
            if ((sideFact * (pos-tempLim)) > 0) {
                tempLim = pos;
            }
            /* even though tempLim moved there is still room on opposing
             * arms for new recombinations, so no need to check armFinished for
             * the opposing arms
             * Since no new armFinished, also no need to recalculate
             * crossoverFinished
             */
        }  
        //update exchangeLim:
        exchangeLim[side] = tempLim;
    } //testRecomb2

    private void testRecomb13(int a) {
        //NOTE: must be called after generateRecombination; 
        //there side, sideFact, pos, mindist are calculated
        //first: check if rejected because of opposing recombination
        //(possibly with interference):
        if (!Double.isNaN(mindist) && rejectDist(mindist)) {
            // conflict with opposing recombination
            armFinished[a] = true;
        }
        else { //test if recombination beyond opposite chromosome end:
            armFinished[a] = (sideFact * (pos-chrom.getSidePos(1-side))) >= 0.0;
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
    } //testRecomb13

    private void makeRecombination(int a) {
        assert(!armFinished[a]);
        recombCount++;
        //insert the recombination into both slots:
        slots[chromatid[0]].insertSegment(pos, chromatid[1]);
        slots[chromatid[1]].insertSegment(pos, chromatid[0]);
        //update lastpos:
        lastpos[a] = pos;
    } //makeRecombination
        

    /**
     * doCrossingOver performs crossing-over for a cross-type Quadrivalent.
     * First we define 8 new HaploStructs ("slots") to represent the chromatids,
     * each with just one segment of one of four pseudo-founderalleles.
     * The four chromosomes are each composed of 2 chromatids ("slots"),
     * 0+1, 2+3, 4+5 and 6+7.
     * Picture the Quadrivalent as a cross, with a left, up, right and down arm.
     * In each arm two of the chromosome ends are paired as follows:
     * left arm (0): head of slots 0,1,6,7
     * up arm (1): tail of slots 0,1,2,3
     * right arm (2): head of slots 2,3,4,5
     * down arm (3): tail of slots 4,5,6,7
     * (in the following we say that the up and down arms are the "opposing
     * arms" of the left and right arms, and v.v.)
     * We define 3 doCrossingOver methods for a cross-type Quadrivalent,
     * which are selected by popdata.quadriMethod = 1/2/3:
     * 1. (implemented in testRecomb13):
     * Recombinations are generated from the ends of each arm, until they
     * conflict with an existing recombinations in an opposite arm or are located
     * beyond the opposite end. When a conflicting recombination has occurred
     * no new recombinations are generated for that arm, but in the other arms
     * the process continues, until in all arms a conflicting recombination has 
     * occurred. This would mean that a conflict affects only the arm from which
     * the next recombination was attempted but not the arm in which the conflicting
     * recombination is located.
     * 2. (implemented in testRecomb2):
     * Recombinations are generated from the ends of each arm until the first
     * conflict occurs. A conflict can occur because:
     * 2a. the new recombination is beyond the opposite end of the chromosome
     * (this is only possible if no recombinations have occurred in the two opposing 
     * arms, else situation 2b applies). - "reason 1" in testRecomb2
     * In this case in effect there are two bivalents, and for one of these
     * bivalents the crossing-over process has just finished.
     * 2b. the new recombination is beyond the nearest recombination in one of the
     * two opposing arms. In that case the chromosome exchange point is now fixed
     * at that recombination position and no new recombinations can be generated
     * in the two arms involved in the conflict. In the other two arms further
     * recombinations are generated until they conflict with the exchange 
     * point. - "reason 2" in testRecomb2
     * 2c. the new recombination is generated before the nearest recombination
     * on the two opposing arm, but close to it and due to interference it is
     * conflicting (this can only occur when popdata.chiasmaInterference is true). 
     * In that case no new recombinations are generated for this arm, but the failed
     * recombination defines the new start point of the chromosome exchange interval. 
     * Also no new recombinations are generated in the arm with the interfering
     * recombination - "reason 3" in testRecomb2
     * 2d. the new recombination is generated beyond the position of a previous
     * conflicting recombination as discussed in 2c, but not conflicting with
     * an actual recombination. Then the start of the exchange interval is
     * moved to that position, which is also the end of the interval -
     * again "reason 2" in testRecomb2
     * The reasoning behind method 2 (a, b, c, d) is that in order to generate
     * a recombination, even if it conflicts, the pairing of the chromosomes in
     * that arm must already have been completed up to that point or up to the 
     * chromosome exchange point.
     * 3. (implemented in testRecomb13):
     * Recombinations are generated from the ends of each arm until the first
     * conflict occurs. When this happens the crossing-over process is
     * ended for the complete quadrivalent.
     * This seems to be too stringent: the amount of recombination at the
     * ends of the arms is about half that in a bivalent.
     * 
     * QUESTION: should the chromosome exchange point cause interference like
     * a recombination? 
     * ANSWER: that depends on whether the interference is mainly caused by the
     * mechanical bending stress (then: yes) or mainly by the break/repair
     * process (then: no). But a practical consideration is that it is often
     * not precisely known where the exchange point is (based on genetics it
     * can only be located to between the two nearest opposing recombination 
     * points), so interference cannot be modelled.
     * Based on simulations with quadrimethod 2 (the default) and without 
     * chiasma interference we see no difference in recombination in small
     * intervals near head, center or tail of the chromosome, but in 200 CM
     * and 400 cM chromosomes we see more recombination in long segments
     * (from about 10% of total chromosome) around the center than at the ends.
     * Perhaps this reflects the more pronounced unimodal distribution of
     * the chromosome exchange point around the chromosome center in the
     * longer chromosomes?
     * 
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
            double recombDist = RECOMBDIST;
            if (!popdata.allowNoRecomb) {
                recombDist = Tools.calcRecombDist(chrom.getLength());
            }
            
            slots = new HaploStruct[8];
            for (int s=0; s<8; s++) {
                slots[s] = new HaploStruct(chrom, -1); //the first segment will be ignored
            }
            for (int a=0; a<4; a++) {
                lastpos[a] = Double.NaN;
                armFinished[a] = false;
            }
            crossoverFinished = false;
            recombCount = 0;
            exchangeLim[0] = chrom.getHeadPos();
            exchangeLim[1] = chrom.getTailPos();
            if (popdata.quadriEachArm) {
                int[] aix; //the sequence of the arms in each cycle
                do {
                    aix = tools.shuffle0123();
                    for (int a=0; a<4; a++) {
                        if (!armFinished[aix[a]]) {
                            generateRecombination(aix[a],recombDist);
                            if (popdata.quadriMethod==2) testRecomb2(aix[a]);
                            else testRecomb13(aix[a]);
                            if (!armFinished[aix[a]]) makeRecombination(aix[a]);
                        } // !armFinished
                    } //for a
                } while (!crossoverFinished); 
            }
            else { //not popdata.quadriEachArm
                //(same as above without aix and the for a - loop)
                do {
                    int a = rand.nextInt(4); //which arm for the next recombination
                        if (!armFinished[a]) {
                            generateRecombination(a, recombDist);
                            if (popdata.quadriMethod==2) testRecomb2(a);
                            else testRecomb13(a);
                            if (!armFinished[a]) makeRecombination(a);
                        } // !armFinished
                } while (!crossoverFinished);
            }
        } while (!(popdata.allowNoRecomb || allChromRecomb(slots)));
        
        if (popdata.testMode) {
            TetraploidChromosome tchrom = (TetraploidChromosome) chrom;
            tchrom.crossQuadrivalentCount++;
            tchrom.quadRecombSum += recombCount;
            tchrom.quadRecombSS += recombCount*recombCount;
            //statistics on midpoint and length of exchange interval:
            if (recombCount==0 ||
                (exchangeLim[0]==chrom.getHeadPos() && 
                 exchangeLim[1]==chrom.getTailPos())) {
                tchrom.noExchangeLimCount++;
            } else if (exchangeLim[0]==chrom.getHeadPos() || 
                       exchangeLim[1]==chrom.getTailPos()) {
                tchrom.oneExchangeLimCount++;
            }
            else {
                tchrom.twoExchangeLimCount++;
                double x = (exchangeLim[0]+exchangeLim[1]) / 2.0; //midpoint of exchange interval
                tchrom.exchangeMidSum += x;
                tchrom.exchangeMidSS += x*x;
                int i = (int)(Math.floor((x-chrom.getHeadPos())*TetraploidChromosome.freqTableLength 
                        / chrom.getLength()));
                tchrom.exchangeMidFreq[i]++;
                x = exchangeLim[1]-exchangeLim[0]; //length of exchange interval
                tchrom.exchangeLengthSum += x;
                tchrom.exchangeLengthSS += x*x;
                i = (int)(Math.floor((x-chrom.getHeadPos())*TetraploidChromosome.freqTableLength 
                        / chrom.getLength()));
                tchrom.exchangeLengthFreq[i]++;
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
            //Note that exchangeLim has been actualized in testRecomb
            boolean centromeresInHeadArm = chrom.getCentromerePos() <= exchangeLim[1]; 
            //adjacent configuration:
            if (centromeresInHeadArm) {
                if (rand.nextBoolean()) 
                     gam = new int[][] { {0,1}, {0,1}, {2,3}, {2,3} };
                else gam = new int[][] { {2,3}, {2,3}, {0,1}, {0,1} };
            }
            else { //centromeres in tail arm
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
