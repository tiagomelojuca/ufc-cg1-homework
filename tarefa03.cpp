// ------------------------------------------------------------------------------------------------
// COMPUTACAO GRAFICA I - TAREFA 03
// ------------------------------------------------------------------------------------------------
// Esfera e planos com sombra
// ------------------------------------------------------------------------------------------------
// Refaca a Tarefa 02 com a inclusao de dois planos. Siga a seguinte especificacao:
// 1) Janela de 60cm x 60 cm (H_J = 60, W_J = 60) 
// 2) Window = Canvas de  500 x 500 pixels (H_C = nLin = 500, W_C = nCol = 500)
// 3) Coordenada z da Janela, z = -d = -30cm
// 4) Esfera com raio R = 40 cm e centro C = (0, 0, - 100cm). 
// 5) Reflectividade da esfera K_d = K_e = K_a = (0.7, 0.2, 0.2),  Shininess = m = 10
// 6) Plano do chao:  Ponto conhecido do plano P_pi = (0, - R, 0),
//                    vetor unitario normal ao plano, n_bar = (0, 1, 0)
// 7) Reflectividade do plano do chao K_d = K_a = (0.2, 0.7, 0.2),
//                                    K_e = (0.0, 0.0, 0.0),  Shininess = m = 1
// 8) Plano de Fundo:  Ponto conhecido do plano P_pi = (0, 0, -200cm),
//                     vetor unitario normal ao plano, n_bar = (0, 0, 1)
// 9) Reflectividade do plano do chao K_d = K_a = (0.3, 0.3, 0.7),
//                                    K_e = (0.0, 0.0, 0.0),  Shininess = m = 1
// 10) Fonte pontual: Intensidade da fonte  I_F = (0.7, 0.7, 0.7),
//                    posicao da fonte P_F = (0, 60cm, -30cm)
// 11) Luz ambiente: Intensidade I_A = (0.3, 0.3, 0.3)

// Obs: 
// 1) Lembre-se que uma intersecao so eh valida se t_i > 0.  
// 2) A cena eh composta de uma esfera, dois planos, uma fonte de luz pontual
//    e uma fonte de luz ambiente. 
// 3) Quando um raio for lancado atraves da janela, ele pode intersectar os tres objetos.
//    Assim o ponto que intersecao visto pelo olho do observador eh aquele
//    que tiver o menor t_i positivo.
// 4) Quando o ponto de intersecao visto estiver no plano do chao ou no plano de fundo,
//    nao calcule as contribuicoes difusa e especular sem antes verificar se
//    o raio P(s) = P_i + s*vetor_l  esta sendo obstruido plea esfera.
// ------------------------------------------------------------------------------------------------
