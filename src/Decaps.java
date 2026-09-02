public class Decaps {

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
                "a.properties"
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

        SecureRandom random = new SecureRandom();


        Element g1 = G0.newRandomElement().getImmutable();
        Element g3 = G0.newRandomElement().getImmutable();
        Element g2 = G1.newRandomElement().getImmutable();


        Element messageM = GT.newRandomElement().getImmutable();
        byte[] gid = new byte[32];
        random.nextBytes(gid);

        int n = 100;
        int repeatCount = 50;


        for (int thetaSize = 3; thetaSize <= 30; thetaSize += 3) {


            int[] Theta = new int[thetaSize];
            for (int sig = 1; sig <= thetaSize; sig++) {
                Theta[sig - 1] = sig;
            }

            Element[] alpha = new Element[thetaSize];
            Element[] x = new Element[thetaSize];
            Element[] A = new Element[thetaSize];
            Element[] B = new Element[thetaSize];
            Element[] C = new Element[thetaSize];


            for (int sig = 1; sig <= thetaSize; sig++) {
                alpha[sig - 1] = Zp.newRandomElement().getImmutable();
            }
            for (int sig = 1; sig <= thetaSize; sig++) {
                x[sig - 1] = Zp.newRandomElement().getImmutable();
            }


            Element pairingBase = pairing
                    .pairing(g1, g2)
                    .getImmutable();


            for (int sig = 1; sig <= thetaSize; sig++) {
                A[sig - 1] = pairingBase
                        .powZn(alpha[sig - 1])
                        .getImmutable();
            }


            for (int sig = 1; sig <= thetaSize; sig++) {
                B[sig - 1] = g2
                        .powZn(x[sig - 1])
                        .getImmutable();
            }


            for (int sig = 1; sig <= thetaSize; sig++) {
                C[sig - 1] = g1
                        .powZn(alpha[sig - 1])
                        .getImmutable();
            }


            for (int l = 10; l <= 100; l += 10) {


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

                Element[][] tldC_i_sig_2 = new Element[l][thetaSize];
                Element[] tldC_i_3 = new Element[l];
                Element[][] tldC_i_sig_4 = new Element[l][thetaSize];
                Element[] tldC_i_5 = new Element[l];



                for (int i = 1; i <= l; i++) {
                    Element e_w_i = pairingBase
                            .powZn(w[i - 1])
                            .getImmutable();

                    for (int sig = 1; sig <= thetaSize; sig++) {
                        Element A_sig_t_i = A[sig - 1]
                                .powZn(t[i - 1])
                                .getImmutable();
                        tldC_i_sig_2[i - 1][sig - 1] = e_w_i
                                .duplicate()
                                .mul(A_sig_t_i)
                                .getImmutable();
                    }
                }


                for (int i = 1; i <= l; i++) {
                    tldC_i_3[i - 1] = g2
                            .powZn(t[i - 1])
                            .getImmutable();
                }


                for (int i = 1; i <= l; i++) {
                    for (int sig = 1; sig <= thetaSize; sig++) {
                        tldC_i_sig_4[i - 1][sig - 1] = B[sig - 1]
                                .powZn(t[i - 1])
                                .getImmutable();
                    }
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
                byte[][] delta = new byte[l][32];
                for (int i = 1; i <= l; i++) {
                    pi[i - 1] = random.nextLong();
                    random.nextBytes(delta[i - 1]);
                }


                Element C1 = tldC1;
                Element[][] C_i_sig_2_1 = new Element[l][thetaSize];
                Element[] C_i_2_2 = new Element[l];
                Element[] C_i_3 = new Element[l];
                Element[][] C_i_sig_4 = new Element[l][thetaSize];
                Element[] C_i_5_1 = new Element[l];
                Element[] C_i_5_2 = new Element[l];
                Element[] C_i_5_3 = new Element[l];


                for (int i = 1; i <= l; i++) {
                    for (int sig = 1; sig <= thetaSize; sig++) {
                        C_i_sig_2_1[i - 1][sig - 1] =
                                tldC_i_sig_2[i - 1][sig - 1];
                    }
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
                    for (int sig = 1; sig <= thetaSize; sig++) {
                        C_i_sig_4[i - 1][sig - 1] =
                                tldC_i_sig_4[i - 1][sig - 1];
                    }
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
                    byte[] piBytes = ByteBuffer
                            .allocate(Long.BYTES)
                            .putLong(pi[i - 1])
                            .array();
                    Element H_pi_i = hashFromBytesToG0(piBytes);
                    Element negative_t_i = t[i - 1]
                            .duplicate()
                            .negate()
                            .getImmutable();
                    C_i_5_3[i - 1] = H_pi_i
                            .powZn(negative_t_i)
                            .getImmutable();
                }


                Element[] r = new Element[thetaSize];
                Element[] atk_sig_1 = new Element[thetaSize];
                Element[] atk_sig_2 = new Element[thetaSize];
                Element[][] atk_i_sig_3 = new Element[l][thetaSize];


                for (int sig = 1; sig <= thetaSize; sig++) {
                    r[sig - 1] = Zp.newRandomElement().getImmutable();
                }


                for (int sig = 1; sig <= thetaSize; sig++) {
                    Element H_gid = hashFromBytesToG0(gid);
                    atk_sig_1[sig - 1] = C[sig - 1]
                            .duplicate()
                            .mul(H_gid.powZn(x[sig - 1]))
                            .getImmutable();
                }


                for (int sig = 1; sig <= thetaSize; sig++) {
                    atk_sig_2[sig - 1] = B[sig - 1]
                            .powZn(r[sig - 1])
                            .getImmutable();
                }


                for (int i = 1; i <= l; i++) {
                    for (int sig = 1; sig <= thetaSize; sig++) {
                        Element H_delta_i = hashFromBytesToG0(delta[i - 1]);
                        atk_i_sig_3[i - 1][sig - 1] = H_delta_i
                                .powZn(r[sig - 1])
                                .getImmutable();
                    }
                }

                Element[] tau = new Element[thetaSize];
                Element[] cek_sig_1 = new Element[thetaSize];
                Element[] cek_sig_2 = new Element[thetaSize];


                for (int sig = 1; sig <= thetaSize; sig++) {
                    tau[sig - 1] = Zp.newRandomElement().getImmutable();
                }


                for (int sig = 1; sig <= thetaSize; sig++) {
                    cek_sig_1[sig - 1] = atk_sig_1[sig - 1]
                            .powZn(tau[sig - 1])
                            .getImmutable();
                }


                for (int sig = 1; sig <= thetaSize; sig++) {
                    cek_sig_2[sig - 1] = g2
                            .powZn(tau[sig - 1])
                            .getImmutable();
                }


                Element[] T = new Element[l];


                for (int i = 1; i <= l; i++) {
                    Element e_g1_g2_C_i_2_2 = pairingBase
                            .powZn(C_i_2_2[i - 1])
                            .getImmutable();

                    Element inverseAttributePairing = pairing
                            .pairing(
                                    atk_i_sig_3[i - 1][0],
                                    C_i_sig_4[i - 1][0]
                            )
                            .invert()
                            .getImmutable();

                    Element combinedC_i_5 = C_i_5_1[i - 1]
                            .duplicate()
                            .mul(g3.powZn(C_i_5_2[i - 1]))
                            .mul(C_i_5_3[i - 1])
                            .getImmutable();

                    Element keyPairing = pairing
                            .pairing(combinedC_i_5, atk_sig_2[0])
                            .getImmutable();

                    T[i - 1] = C_i_sig_2_1[i - 1][0]
                            .duplicate()
                            .mul(e_g1_g2_C_i_2_2)
                            .mul(inverseAttributePairing)
                            .mul(keyPairing)
                            .getImmutable();
                }


                Element[] d = new Element[l];
                for (int i = 1; i <= l; i++) {
                    d[i - 1] = Zp.newRandomElement().getImmutable();
                }

                double totalTimeMs = 0.0;
                long totalMBytes = 0;


                for (int repeat = 0; repeat < repeatCount; repeat++) {

                    long startTime = System.nanoTime();


                    Element Temp = GT.newOneElement();
                    for (int i = 1; i <= l; i++) {
                        Element H_gid = hashFromBytesToG0(gid);
                        Element gidPairing = pairing
                                .pairing(H_gid, C_i_sig_4[i - 1][0])
                                .getImmutable();

                        Element factor = T[i - 1]
                                .duplicate()
                                .mul(gidPairing)
                                .getImmutable();

                        Element negative_d_i = d[i - 1]
                                .duplicate()
                                .negate()
                                .getImmutable();

                        Element weightedFactor = factor
                                .powZn(negative_d_i)
                                .getImmutable();

                        Temp.mul(weightedFactor);
                    }
                    Temp = Temp.getImmutable();


                    Element recoveredM = C1
                            .duplicate()
                            .mul(Temp)
                            .getImmutable();

                    long endTime = System.nanoTime();
                    totalTimeMs += (endTime - startTime) / 1_000_000.0;


                    totalMBytes += recoveredM.toBytes().length;
                }

                double averageTimeMs = totalTimeMs / repeatCount;
                double averageMBytes = (double) totalMBytes / repeatCount;

                System.out.printf("Theta=%d, l=%d%n", thetaSize, l);
                System.out.printf(
                        "Decaps computation time: %.3f ms%n",
                        averageTimeMs
                );
                System.out.printf(
                        "M storage: %.3f Bytes%n",
                        averageMBytes
                );
                System.out.println("--------------------------------------------------");


                if (C_i_3[l - 1] == null
                        || atk_sig_1[thetaSize - 1] == null
                        || cek_sig_1[thetaSize - 1] == null
                        || cek_sig_2[thetaSize - 1] == null) {
                    throw new IllegalStateException("Input preparation failed.");
                }
            }

            if (Theta.length != thetaSize
                    || A.length != thetaSize
                    || B.length != thetaSize
                    || C.length != thetaSize) {
                throw new IllegalStateException("Authority initialization failed.");
            }
        }
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
