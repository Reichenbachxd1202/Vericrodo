public class UtkGen {

    static Pairing pairing;
    static Field<Element> G0, G1, GT, Zp;


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

                double totalType1TimeMs = 0.0;
                double totalType2TimeMs = 0.0;
                double totalType3TimeMs = 0.0;
                long totalTok1Bytes = 0;
                long totalTok2Bytes = 0;
                long totalTok3Bytes = 0;

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
                        Element e_w_i = pairingBase
                                .powZn(w[i - 1])
                                .getImmutable();
                        Element A_1_t_i = A_1
                                .powZn(t[i - 1])
                                .getImmutable();
                        tldC_i_1_2[i - 1] = e_w_i
                                .duplicate()
                                .mul(A_1_t_i)
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


                    Element[] vecW = multiplyMatrixByVector(shareMatrix, vecV);


                    Element[] vecEta = new Element[n];
                    vecEta[0] = Zp.newZeroElement().getImmutable();
                    for (int j = 2; j <= n; j++) {
                        vecEta[j - 1] = Zp.newRandomElement().getImmutable();
                    }


                    Element[] vecXi = multiplyMatrixByVector(shareMatrix, vecEta);


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
                        byte[] piBytes = longToBytes(pi[i - 1]);
                        Element H_pi_i = hashFromBytesToG0(piBytes);
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




                    Element[] uk1_i_1 = new Element[l];
                    Element[] uk1_i_2 = new Element[l];
                    Element[] uk1_i_3 = new Element[l];
                    Element[] uk1_i_4 = new Element[l];
                    Element[] uk1_i_5 = new Element[l];
                    Element[] uk1_i_6 = new Element[l];

                    long type1StartTime = System.nanoTime();


                    for (int i = 1; i <= l; i++) {
                        uk1_i_1[i - 1] = subtractProduct(
                                wpr[i - 1],
                                beta[i - 1],
                                w[i - 1]
                        );
                    }


                    for (int i = 1; i <= l; i++) {
                        uk1_i_2[i - 1] = subtractProduct(
                                tpr[i - 1],
                                beta[i - 1],
                                t[i - 1]
                        );
                    }


                    for (int i = 1; i <= l; i++) {
                        uk1_i_3[i - 1] = subtractProduct(
                                vecWpr[i - 1],
                                beta[i - 1],
                                vecW[i - 1]
                        );
                    }


                    for (int i = 1; i <= l; i++) {
                        uk1_i_4[i - 1] = subtractProduct(
                                xipr[i - 1],
                                beta[i - 1],
                                xi[i - 1]
                        );
                    }


                    for (int i = 1; i <= l; i++) {
                        uk1_i_5[i - 1] = subtractProduct(
                                vecXipr[i - 1],
                                beta[i - 1],
                                vecXi[i - 1]
                        );
                    }


                    for (int i = 1; i <= l; i++) {
                        Element H_pipr_i = hashFromBytesToG0(
                                longToBytes(pipr[i - 1])
                        );
                        Element negative_uk1_i_2 = uk1_i_2[i - 1]
                                .duplicate()
                                .negate()
                                .getImmutable();
                        uk1_i_6[i - 1] = H_pipr_i
                                .powZn(negative_uk1_i_2)
                                .getImmutable();
                    }

                    long type1EndTime = System.nanoTime();
                    totalType1TimeMs += (type1EndTime - type1StartTime)
                            / 1_000_000.0;

                    long tok1Bytes = 0;
                    for (int i = 1; i <= l; i++) {
                        tok1Bytes += uk1_i_1[i - 1].toBytes().length;
                        tok1Bytes += uk1_i_2[i - 1].toBytes().length;
                        tok1Bytes += uk1_i_3[i - 1].toBytes().length;
                        tok1Bytes += uk1_i_4[i - 1].toBytes().length;
                        tok1Bytes += uk1_i_5[i - 1].toBytes().length;
                        tok1Bytes += uk1_i_6[i - 1].toBytes().length;
                    }
                    totalTok1Bytes += tok1Bytes;




                    Element[] uk2_i_1 = new Element[l];
                    Element[] uk2_i_2 = new Element[l];
                    Element[] uk2_i_3 = new Element[l];
                    Element[] uk2_i_4 = new Element[l];
                    Element[] uk2_i_5 = new Element[l];

                    long type2StartTime = System.nanoTime();


                    for (int i = 1; i <= l; i++) {
                        uk2_i_1[i - 1] = subtractProduct(
                                wpr[i - 1],
                                beta[i - 1],
                                w[i - 1]
                        );
                    }


                    for (int i = 1; i <= l; i++) {
                        uk2_i_2[i - 1] = subtractProduct(
                                t[i - 1],
                                beta[i - 1],
                                t[i - 1]
                        );
                    }


                    for (int i = 1; i <= l; i++) {
                        uk2_i_3[i - 1] = subtractProduct(
                                vecWpr[i - 1],
                                beta[i - 1],
                                vecW[i - 1]
                        );
                    }


                    for (int i = 1; i <= l; i++) {
                        uk2_i_4[i - 1] = subtractProduct(
                                xipr[i - 1],
                                beta[i - 1],
                                xi[i - 1]
                        );
                    }


                    for (int i = 1; i <= l; i++) {
                        uk2_i_5[i - 1] = subtractProduct(
                                vecXipr[i - 1],
                                beta[i - 1],
                                vecXi[i - 1]
                        );
                    }

                    long type2EndTime = System.nanoTime();
                    totalType2TimeMs += (type2EndTime - type2StartTime)
                            / 1_000_000.0;

                    long tok2Bytes = 0;
                    for (int i = 1; i <= l; i++) {
                        tok2Bytes += uk2_i_1[i - 1].toBytes().length;
                        tok2Bytes += uk2_i_2[i - 1].toBytes().length;
                        tok2Bytes += uk2_i_3[i - 1].toBytes().length;
                        tok2Bytes += uk2_i_4[i - 1].toBytes().length;
                        tok2Bytes += uk2_i_5[i - 1].toBytes().length;
                    }
                    totalTok2Bytes += tok2Bytes;




                    Element[] uk3_i_1 = new Element[l];
                    Element[] uk3_i_2 = new Element[l];
                    Element[] uk3_i_3 = new Element[l];
                    Element[] uk3_i_4 = new Element[l];
                    Element[] uk3_i_5 = new Element[l];
                    Element[] uk3_i_6 = new Element[l];
                    Element[] uk3_i_7 = new Element[l];

                    long type3StartTime = System.nanoTime();


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

                    long type3EndTime = System.nanoTime();
                    totalType3TimeMs += (type3EndTime - type3StartTime)
                            / 1_000_000.0;

                    long tok3Bytes = 0;
                    for (int i = 1; i <= l; i++) {
                        tok3Bytes += uk3_i_1[i - 1].toBytes().length;
                        tok3Bytes += uk3_i_2[i - 1].toBytes().length;
                        tok3Bytes += uk3_i_3[i - 1].toBytes().length;
                        tok3Bytes += uk3_i_4[i - 1].toBytes().length;
                        tok3Bytes += uk3_i_5[i - 1].toBytes().length;
                        tok3Bytes += uk3_i_6[i - 1].toBytes().length;
                        tok3Bytes += uk3_i_7[i - 1].toBytes().length;
                    }
                    totalTok3Bytes += tok3Bytes;


                    if (C1 == null
                            || C_i_1_2_1[l - 1] == null
                            || C_i_2_2[l - 1] == null
                            || C_i_3[l - 1] == null
                            || C_i_1_4[l - 1] == null
                            || C_i_5_1[l - 1] == null
                            || C_i_5_2[l - 1] == null
                            || C_i_5_3[l - 1] == null) {
                        throw new IllegalStateException(
                                "OnEncaps preparation failed."
                        );
                    }
                }

                double averageType1TimeMs = totalType1TimeMs / repeatCount;
                double averageType2TimeMs = totalType2TimeMs / repeatCount;
                double averageType3TimeMs = totalType3TimeMs / repeatCount;
                double averageTok1KBytes = (double) totalTok1Bytes
                        / repeatCount
                        / 1024.0;
                double averageTok2KBytes = (double) totalTok2Bytes
                        / repeatCount
                        / 1024.0;
                double averageTok3KBytes = (double) totalTok3Bytes
                        / repeatCount
                        / 1024.0;

                System.out.printf("Theta=%d, l=%d, sig=1%n", thetaSize, l);
                System.out.printf(
                        "Type 1 Computation Time: %.3f ms%n",
                        averageType1TimeMs
                );
                System.out.printf(
                        "Tok1 storage: %.3f KBytes%n",
                        averageTok1KBytes
                );
                System.out.printf(
                        "Type 2 Computation Time: %.3f ms%n",
                        averageType2TimeMs
                );
                System.out.printf(
                        "Tok2 storage: %.3f KBytes%n",
                        averageTok2KBytes
                );
                System.out.printf(
                        "Type 3 Computation Time: %.3f ms%n",
                        averageType3TimeMs
                );
                System.out.printf(
                        "Tok3 storage: %.3f KBytes%n",
                        averageTok3KBytes
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

    private static Element subtractProduct(
            Element minuend,
            Element factor1,
            Element factor2
    ) {
        Element product = factor1
                .duplicate()
                .mul(factor2)
                .getImmutable();
        return minuend
                .duplicate()
                .sub(product)
                .getImmutable();
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
