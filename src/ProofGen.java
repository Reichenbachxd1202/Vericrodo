public class ProofGen {

    static Pairing pairing;
    static Field<Element> G0, G1, GT, Zp;
    static volatile Object benchmarkSink;


    public static Element hashFromBytesToG0(byte[] bytes) {
        return G0.newElement()
                .setFromHash(bytes, 0, bytes.length)
                .getImmutable();
    }


    public static int hashFromUToTheta(Element element, int[] Theta) {
        if (Theta == null || Theta.length == 0) {
            throw new IllegalArgumentException("Theta must not be empty.");
        }

        try {
            MessageDigest digest = MessageDigest.getInstance("SHA-256");
            byte[] hash = digest.digest(element.toBytes());
            int index = new BigInteger(1, hash)
                    .mod(BigInteger.valueOf(Theta.length))
                    .intValue();
            return Theta[index];
        } catch (NoSuchAlgorithmException e) {
            throw new RuntimeException("SHA-256 is not available.", e);
        }
    }

    public static void main(String[] args) {


        List<String> curveList = Arrays.asList(
                "g149.properties",
                "f.properties"
        );

        for (String curveFile : curveList) {
            System.out.printf("Curve: %s%n", curveFile);
            runForCurve(curveFile);
            System.out.println("==================================================");
        }
    }

    private static void runForCurve(String curveFile) {

        pairing = PairingFactory.getPairing(curveFile);


        G0 = pairing.getG1();
        G1 = pairing.getG2();
        GT = pairing.getGT();
        Zp = pairing.getZr();

        Random random = new Random();


        Element g1 = G0.newRandomElement().getImmutable();
        Element g3 = G0.newRandomElement().getImmutable();
        Element g2 = G1.newRandomElement().getImmutable();


        Element messageM = GT.newRandomElement().getImmutable();

        int repeatCount = 20;
        int n = 100;


        for (int thetaSize = 3; thetaSize <= 30; thetaSize += 3) {


            int[] Theta = new int[thetaSize];
            for (int authorityId = 1;
                    authorityId <= thetaSize;
                    authorityId++) {
                Theta[authorityId - 1] = authorityId;
            }


            Element alpha_1 = Zp.newRandomElement().getImmutable();
            Element x_1 = Zp.newRandomElement().getImmutable();


            Element pairingBase = pairing
                    .pairing(g1, g2)
                    .getImmutable();


            Element A_1 = pairingBase
                    .powZn(alpha_1)
                    .getImmutable();


            Element B_1 = g2
                    .powZn(x_1)
                    .getImmutable();


            Element authorityC_1 = g1
                    .powZn(alpha_1)
                    .getImmutable();


            for (int l = 10; l <= 100; l += 10) {

                double totalProofGenTimeMs = 0.0;
                long totalProofBytes = 0;

                for (int repeat = 0; repeat < repeatCount; repeat++) {




                    Element v = Zp.newRandomElement().getImmutable();
                    Element[] w = new Element[l];
                    Element[] xi = new Element[l];
                    Element[] t = new Element[l];

                    for (int i = 1; i <= l; i++) {
                        w[i - 1] = Zp.newRandomElement().getImmutable();
                    }
                    for (int i = 1; i <= l; i++) {
                        xi[i - 1] = Zp.newRandomElement().getImmutable();
                    }
                    for (int i = 1; i <= l; i++) {
                        t[i - 1] = Zp.newRandomElement().getImmutable();
                    }


                    Element tldC1 = messageM
                            .duplicate()
                            .mul(pairingBase.powZn(v))
                            .getImmutable();

                    Element[] tldC_i_1_2 = new Element[l];
                    Element[] tldC_i_3 = new Element[l];
                    Element[] tldC_i_1_4 = new Element[l];
                    Element[] tldC_i_5 = new Element[l];


                    for (int i = 1; i <= l; i++) {
                        Element firstTerm = pairingBase
                                .powZn(w[i - 1])
                                .getImmutable();
                        Element secondTerm = A_1
                                .powZn(t[i - 1])
                                .getImmutable();
                        tldC_i_1_2[i - 1] = firstTerm
                                .duplicate()
                                .mul(secondTerm)
                                .getImmutable();
                    }


                    for (int i = 1; i <= l; i++) {
                        tldC_i_3[i - 1] = g2
                                .powZn(t[i - 1])
                                .getImmutable();
                    }


                    for (int i = 1; i <= l; i++) {
                        tldC_i_1_4[i - 1] = B_1
                                .powZn(t[i - 1])
                                .getImmutable();
                    }


                    for (int i = 1; i <= l; i++) {
                        tldC_i_5[i - 1] = g3
                                .powZn(xi[i - 1])
                                .getImmutable();
                    }


                    Element[][] shareMatrix = new Element[l][n];
                    for (int i = 1; i <= l; i++) {
                        for (int j = 1; j <= n; j++) {
                            shareMatrix[i - 1][j - 1] = Zp
                                    .newRandomElement()
                                    .getImmutable();
                        }
                    }


                    Element[] vecV = new Element[n];
                    vecV[0] = v;
                    for (int j = 2; j <= n; j++) {
                        vecV[j - 1] = Zp.newRandomElement().getImmutable();
                    }


                    Element[] vecW = multiplyMatrixByVector(
                            shareMatrix, vecV
                    );


                    Element[] vecEta = new Element[n];
                    vecEta[0] = Zp.newZeroElement().getImmutable();
                    for (int j = 2; j <= n; j++) {
                        vecEta[j - 1] = Zp.newRandomElement().getImmutable();
                    }


                    Element[] vecXi = multiplyMatrixByVector(
                            shareMatrix, vecEta
                    );


                    long[] pi = new long[l];
                    for (int i = 1; i <= l; i++) {
                        pi[i - 1] = random.nextLong();
                    }


                    Element C1 = tldC1;
                    Element[] C_i_1_2_1 = new Element[l];
                    Element[] C_i_2_2 = new Element[l];
                    Element[] C_i_3 = new Element[l];
                    Element[] C_i_1_4 = new Element[l];
                    Element[] C_i_5_1 = new Element[l];
                    Element[] C_i_5_2 = new Element[l];
                    Element[] C_i_5_3 = new Element[l];


                    for (int i = 1; i <= l; i++) {
                        C_i_1_2_1[i - 1] = tldC_i_1_2[i - 1];
                    }


                    for (int i = 1; i <= l; i++) {
                        C_i_2_2[i - 1] = vecW[i - 1]
                                .duplicate()
                                .sub(w[i - 1])
                                .getImmutable();
                    }


                    for (int i = 1; i <= l; i++) {
                        C_i_3[i - 1] = tldC_i_3[i - 1];
                    }


                    for (int i = 1; i <= l; i++) {
                        C_i_1_4[i - 1] = tldC_i_1_4[i - 1];
                    }


                    for (int i = 1; i <= l; i++) {
                        C_i_5_1[i - 1] = tldC_i_5[i - 1];
                    }


                    for (int i = 1; i <= l; i++) {
                        C_i_5_2[i - 1] = vecXi[i - 1]
                                .duplicate()
                                .sub(xi[i - 1])
                                .getImmutable();
                    }


                    for (int i = 1; i <= l; i++) {
                        Element H_pi_i = hashFromBytesToG0(
                                longToBytes(pi[i - 1])
                        );
                        Element negative_t_i = t[i - 1]
                                .duplicate()
                                .negate()
                                .getImmutable();
                        C_i_5_3[i - 1] = H_pi_i
                                .powZn(negative_t_i)
                                .getImmutable();
                    }




                    Element[] beta = new Element[l];
                    Element[] wpr = new Element[l];
                    Element[] tpr = new Element[l];
                    Element[] xipr = new Element[l];
                    Element[] vecWpr = new Element[l];
                    Element[] vecXipr = new Element[l];
                    long[] pipr = new long[l];

                    for (int i = 1; i <= l; i++) {
                        beta[i - 1] = Zp.newRandomElement().getImmutable();
                        wpr[i - 1] = Zp.newRandomElement().getImmutable();
                        tpr[i - 1] = Zp.newRandomElement().getImmutable();
                        xipr[i - 1] = Zp.newRandomElement().getImmutable();
                    }

                    for (int i = 1; i <= l; i++) {
                        vecWpr[i - 1] = Zp.newRandomElement().getImmutable();
                        vecXipr[i - 1] = Zp.newRandomElement().getImmutable();
                    }

                    for (int i = 1; i <= l; i++) {
                        pipr[i - 1] = random.nextLong();
                    }

                    Element[] uk3_i_1 = new Element[l];
                    Element[] uk3_i_2 = new Element[l];
                    Element[] uk3_i_3 = new Element[l];
                    Element[] uk3_i_4 = new Element[l];
                    Element[] uk3_i_5 = new Element[l];
                    Element[] uk3_i_6 = new Element[l];
                    Element[] uk3_i_7 = new Element[l];


                    for (int i = 1; i <= l; i++) {
                        Element firstTerm = pairingBase
                                .powZn(wpr[i - 1])
                                .getImmutable();
                        Element secondTerm = A_1
                                .powZn(tpr[i - 1])
                                .getImmutable();
                        uk3_i_1[i - 1] = firstTerm
                                .duplicate()
                                .mul(secondTerm)
                                .getImmutable();
                    }


                    for (int i = 1; i <= l; i++) {
                        uk3_i_2[i - 1] = vecWpr[i - 1]
                                .duplicate()
                                .sub(wpr[i - 1])
                                .getImmutable();
                    }


                    for (int i = 1; i <= l; i++) {
                        uk3_i_3[i - 1] = g2
                                .powZn(tpr[i - 1])
                                .getImmutable();
                    }


                    for (int i = 1; i <= l; i++) {
                        uk3_i_4[i - 1] = B_1
                                .powZn(tpr[i - 1])
                                .getImmutable();
                    }


                    for (int i = 1; i <= l; i++) {
                        uk3_i_5[i - 1] = g3
                                .powZn(xipr[i - 1])
                                .getImmutable();
                    }


                    for (int i = 1; i <= l; i++) {
                        uk3_i_6[i - 1] = vecXipr[i - 1]
                                .duplicate()
                                .sub(xipr[i - 1])
                                .getImmutable();
                    }


                    for (int i = 1; i <= l; i++) {
                        Element H_pipr_i = hashFromBytesToG0(
                                longToBytes(pipr[i - 1])
                        );
                        Element negative_tpr_i = tpr[i - 1]
                                .duplicate()
                                .negate()
                                .getImmutable();
                        uk3_i_7[i - 1] = H_pipr_i
                                .powZn(negative_tpr_i)
                                .getImmutable();
                    }




                    Element[] Cpr_i_1_2_1 = new Element[l];
                    Element[] Cpr_i_2_2 = new Element[l];
                    Element[] Cpr_i_3 = new Element[l];
                    Element[] Cpr_i_1_4 = new Element[l];
                    Element[] Cpr_i_5_1 = new Element[l];
                    Element[] Cpr_i_5_2 = new Element[l];
                    Element[] Cpr_i_5_3 = new Element[l];

                    for (int i = 1; i <= l; i++) {
                        Cpr_i_1_2_1[i - 1] = uk3_i_1[i - 1];
                    }
                    for (int i = 1; i <= l; i++) {
                        Cpr_i_2_2[i - 1] = uk3_i_2[i - 1];
                    }
                    for (int i = 1; i <= l; i++) {
                        Cpr_i_3[i - 1] = uk3_i_3[i - 1];
                    }
                    for (int i = 1; i <= l; i++) {
                        Cpr_i_1_4[i - 1] = uk3_i_4[i - 1];
                    }
                    for (int i = 1; i <= l; i++) {
                        Cpr_i_5_1[i - 1] = uk3_i_5[i - 1];
                    }
                    for (int i = 1; i <= l; i++) {
                        Cpr_i_5_2[i - 1] = uk3_i_6[i - 1];
                    }
                    for (int i = 1; i <= l; i++) {
                        Cpr_i_5_3[i - 1] = uk3_i_7[i - 1];
                    }


                    Element[] pf = new Element[l];




                    long startTime = System.nanoTime();


                    Element gamma = Zp.newRandomElement().getImmutable();


                    Element g2Gamma = g2
                            .powZn(gamma)
                            .getImmutable();


                    for (int i = 1; i <= l; i++) {
                        Element firstTerm = Cpr_i_1_2_1[i - 1];

                        Element secondTerm = pairingBase
                                .powZn(Cpr_i_2_2[i - 1])
                                .getImmutable();

                        Element H_pipr_i = hashFromBytesToG0(
                                longToBytes(pipr[i - 1])
                        );
                        Element H_pipr_i_gamma = H_pipr_i
                                .powZn(gamma)
                                .getImmutable();
                        Element thirdTerm = pairing
                                .pairing(
                                        H_pipr_i_gamma,
                                        Cpr_i_3[i - 1]
                                )
                                .getImmutable();

                        Element g3ExponentTerm = g3
                                .powZn(Cpr_i_5_2[i - 1])
                                .getImmutable();
                        Element fourthPairingInput = Cpr_i_5_1[i - 1]
                                .duplicate()
                                .mul(g3ExponentTerm)
                                .mul(Cpr_i_5_3[i - 1])
                                .getImmutable();
                        Element fourthTerm = pairing
                                .pairing(fourthPairingInput, g2Gamma)
                                .getImmutable();

                        pf[i - 1] = firstTerm
                                .duplicate()
                                .mul(secondTerm)
                                .mul(thirdTerm)
                                .mul(fourthTerm)
                                .getImmutable();
                    }

                    long endTime = System.nanoTime();
                    totalProofGenTimeMs += (endTime - startTime)
                            / 1_000_000.0;


                    long proofBytes = 0;
                    for (int i = 1; i <= l; i++) {
                        proofBytes += pf[i - 1].toBytes().length;
                        proofBytes += Cpr_i_1_4[i - 1].toBytes().length;
                    }
                    totalProofBytes += proofBytes;


                    benchmarkSink = new Object[] {
                            pf, Cpr_i_1_4, gamma, C1, beta
                    };

                    if (pf[l - 1] == null
                            || Cpr_i_1_4[l - 1] == null
                            || C1 == null
                            || beta[l - 1] == null) {
                        throw new IllegalStateException(
                                "Proof generation failed."
                        );
                    }
                }

                double averageProofGenTimeMs = totalProofGenTimeMs
                        / repeatCount;
                double averageProofKBytes = (double) totalProofBytes
                        / repeatCount
                        / 1024.0;

                System.out.printf(
                        "Theta=%d, l=%d, sig=1%n",
                        thetaSize, l
                );
                System.out.printf(
                        "ProofGen Computation Time: %.3f ms%n",
                        averageProofGenTimeMs
                );
                System.out.printf(
                        "PF storage: %.3f KBytes%n",
                        averageProofKBytes
                );
                System.out.println("--------------------------------------------------");
            }

            if (Theta.length != thetaSize
                    || authorityC_1 == null) {
                throw new IllegalStateException(
                        "Theta or sig=1 initialization failed."
                );
            }
        }
    }

    private static byte[] longToBytes(long value) {
        return ByteBuffer
                .allocate(Long.BYTES)
                .putLong(value)
                .array();
    }

    private static Element[] multiplyMatrixByVector(
            Element[][] matrix,
            Element[] vector
    ) {
        int rows = matrix.length;
        int columns = vector.length;
        Element[] result = new Element[rows];

        for (int i = 0; i < rows; i++) {
            Element sum = Zp.newZeroElement();
            for (int j = 0; j < columns; j++) {
                sum.add(matrix[i][j].duplicate().mul(vector[j]));
            }
            result[i] = sum.getImmutable();
        }

        return result;
    }
}
