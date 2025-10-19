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

#include <cstdint>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <memory>

#include "BitmapPlusPlus.hpp"

// ------------------------------------------------------------------------------------------------

class TVetor4D
{
public:
    TVetor4D() = default;
    TVetor4D(double x, double y, double z, bool pt = true) : _x(x), _y(y), _z(z), _p(pt) {}

    double X() const { return _x; }
    double Y() const { return _y; }
    double Z() const { return _z; }

    void X(double x) { _x = x; }
    void Y(double y) { _y = y; }
    void Z(double z) { _z = z; }

    TVetor4D& operator+=(const TVetor4D& outro)
    {
        _x += outro._x;
        _y += outro._y;
        _z += outro._z;
        _p += outro._p;
        return *this;
    }

    TVetor4D operator+(const TVetor4D& outro) const
    {
        TVetor4D ret = *this;
        ret += outro;
        return ret;
    }

    TVetor4D& operator*=(double k)
    {
        _x *= k;
        _y *= k;
        _z *= k;
        _p *= k;
        return *this;
    }

    TVetor4D operator*(double k) const
    {
        TVetor4D ret = *this;
        ret *= k;
        return ret;
    }

    TVetor4D& operator-=(const TVetor4D& outro)
    {
        TVetor4D copiaOutro = outro;
        copiaOutro *= -1.0;
        *this += copiaOutro;
        return *this;
    }

    TVetor4D operator-(const TVetor4D& outro) const
    {
        TVetor4D ret = *this;
        ret -= outro;
        return ret;
    }

    TVetor4D& operator/=(double k)
    {
        _x /= k;
        _y /= k;
        _z /= k;
        _p /= k;
        return *this;
    }

    TVetor4D operator/(double k) const
    {
        TVetor4D ret = *this;
        ret /= k;
        return ret;
    }

    void Clamp(double min, double max)
    {
        ClampX(min, max);
        ClampY(min, max);
        ClampZ(min, max);
    }
    void ClampX(double min, double max)
    {
        _x = std::clamp(_x, min, max);
    }
    void ClampY(double min, double max)
    {
        _y = std::clamp(_y, min, max);
    }
    void ClampZ(double min, double max)
    {
        _z = std::clamp(_z, min, max);
    }

protected:
    double _x = 0.0;
    double _y = 0.0;
    double _z = 0.0;
    int8_t _p =   1;
};

// ------------------------------------------------------------------------------------------------

class TPonto3D : public TVetor4D
{
public:
    TPonto3D() = default;
    TPonto3D(double x, double y, double z) : TVetor4D(x, y, z) {}
    TPonto3D(const TVetor4D& outro) : TVetor4D(outro.X(), outro.Y(), outro.Z()) {}
};

// ------------------------------------------------------------------------------------------------

class TVetor3D : public TVetor4D
{
public:
    TVetor3D() : TVetor4D(0.0, 0.0, 0.0, false) {}
    TVetor3D(double x, double y, double z) : TVetor4D(x, y, z, false) {}
    TVetor3D(const TVetor4D& outro) : TVetor4D(outro.X(), outro.Y(), outro.Z(), false) {}

    double Dot(const TVetor3D& outro) const
    {
        return _x * outro._x +
               _y * outro._y +
               _z * outro._z;
    }

    TVetor3D Arroba(const TVetor3D& outro) const
    {
        return {
            _x * outro._x,
            _y * outro._y,
            _z * outro._z
        };
    }

    double Norma() const
    {
        return sqrt(Dot(*this));
    }

    TVetor3D Normalizado() const
    {
        TVetor3D ret = *this;
        ret /= Norma();
        return ret;
    }
};

// ------------------------------------------------------------------------------------------------

class TRaio3D
{
public:
    TRaio3D() = default;
    TRaio3D(const TPonto3D& origem, const TVetor3D& direcao) : _origem(origem), _direcao(direcao) {}

    const TPonto3D& Origem() const
    {
        return _origem;
    }
    const TVetor3D& Direcao() const
    {
        return _direcao;
    }
    TPonto3D Ponto(double t) const
    {
        return _origem + _direcao * t;
    }

private:
    TPonto3D _origem;
    TVetor3D _direcao;
};

// ------------------------------------------------------------------------------------------------

class TCor
{
public:
    TCor() = default;
    TCor(uint8_t r, uint8_t g, uint8_t b) : _r(r), _g(g), _b(b) {}

    bool operator==(const TCor& outra) const
    {
        return _r == outra._r && _g == outra._g && _b == outra._b;
    }

    uint8_t R() const { return _r; };
    uint8_t G() const { return _g; };
    uint8_t B() const { return _b; };

    void R(uint8_t r) { _r = r; };
    void G(uint8_t g) { _g = g; };
    void B(uint8_t b) { _b = b; };

private:
    uint8_t _r = 0;
    uint8_t _g = 0;
    uint8_t _b = 0;
};

// ------------------------------------------------------------------------------------------------

class TMaterial
{
public:
    TMaterial() = default;

    double M() const
    {
        return _m;
    }
    void M(double m)
    {
        _m = m;
    }

    double KdR() const
    {
        return _kdR;
    }
    void KdR(double r)
    {
        _kdR = r;
    }
    double KdG() const
    {
        return _kdG;
    }
    void KdG(double g)
    {
        _kdG = g;
    }
    double KdB() const
    {
        return _kdB;
    }
    void KdB(double b)
    {
        _kdB = b;
    }

    double KeR() const
    {
        return _keR;
    }
    void KeR(double r)
    {
        _keR = r;
    }
    double KeG() const
    {
        return _keG;
    }
    void KeG(double g)
    {
        _keG = g;
    }
    double KeB() const
    {
        return _keB;
    }
    void KeB(double b)
    {
        _keB = b;
    }

    double KaR() const
    {
        return _kaR;
    }
    void KaR(double r)
    {
        _kaR = r;
    }
    double KaG() const
    {
        return _kaG;
    }
    void KaG(double g)
    {
        _kaG = g;
    }
    double KaB() const
    {
        return _kaB;
    }
    void KaB(double b)
    {
        _kaB = b;
    }

private:
    double _m = 0.0;

    double _kdR = 0.0;
    double _kdG = 0.0;
    double _kdB = 0.0;

    double _keR = 0.0;
    double _keG = 0.0;
    double _keB = 0.0;

    double _kaR = 0.0;
    double _kaG = 0.0;
    double _kaB = 0.0;
};

// ------------------------------------------------------------------------------------------------

class IFonteLuminosa
{
public:
    IFonteLuminosa() = default;
    virtual ~IFonteLuminosa() = default;

    virtual IFonteLuminosa* Copia() const = 0;
};

// ------------------------------------------------------------------------------------------------

class TFontePontual : public IFonteLuminosa
{
public:
    TFontePontual() = default;
    TFontePontual(const TPonto3D& posicao, const TVetor3D& intensidade)
        : _p(posicao), _i(intensidade) {}

    IFonteLuminosa* Copia() const override
    {
        return new TFontePontual(*this);
    }

    const TPonto3D& Posicao() const
    {
        return _p;
    }
    const TVetor3D& Intensidade() const
    {
        return _i;
    }

private:
    TPonto3D _p;
    TVetor3D _i;
};

// ------------------------------------------------------------------------------------------------

class IEntidade3D
{
public:
    IEntidade3D() = default;
    virtual ~IEntidade3D() = default;

    virtual IEntidade3D* Copia() const = 0;
    virtual std::string Rotulo() const = 0;

    virtual std::vector<double> Intersecoes(const TRaio3D& raio) const = 0;
};

// ------------------------------------------------------------------------------------------------

class TPlano : public IEntidade3D
{
public:
    TPlano() = delete;
    TPlano(const TPonto3D& posicao, const TVetor3D& direcao) : _p(posicao), _n(direcao) {}

    IEntidade3D* Copia() const override
    {
        return new TPlano(*this);
    }

    std::string Rotulo() const override
    {
        return _rotulo;
    }
    void Rotulo(const std::string& rotulo)
    {
        _rotulo = rotulo;
    }

    const TPonto3D& PontoReferencia() const
    {
        return _p;
    }
    const TVetor3D& Normal() const
    {
        return _n;
    }

    const TMaterial& Material() const { return _material; }
    void Material(const TMaterial& material) { _material = material; }

    std::vector<double> Intersecoes(const TRaio3D& raio) const override
    {
        std::vector<double> intersecoes;

        const double dn = raio.Direcao().Dot(_n);
        if (dn != 0.0)
        {
            intersecoes.push_back(-(_n.Dot(raio.Origem() - _p) / dn));
        }

        return intersecoes;
    }

private:
    std::string _rotulo;

    TPonto3D _p;
    TVetor3D _n;

    TMaterial _material;
};

// ------------------------------------------------------------------------------------------------

class TEsfera : public IEntidade3D
{
public:
    TEsfera() = delete;
    TEsfera(const TPonto3D& centro, double raio) : _centro(centro), _raio(raio) {}

    IEntidade3D* Copia() const override
    {
        return new TEsfera(*this);
    }

    std::string Rotulo() const override
    {
        return _rotulo;
    }
    void Rotulo(const std::string& rotulo)
    {
        _rotulo = rotulo;
    }

    const TPonto3D& Centro() const { return _centro; }
    double Raio() const { return _raio; }

    const TMaterial& Material() const { return _material; }
    void Material(const TMaterial& material) { _material = material; }

    std::vector<double> Intersecoes(const TRaio3D& raio) const override
    {
        std::vector<double> intersecoes;

        const TVetor3D w = raio.Origem() - Centro();
        const TVetor3D d = raio.Direcao();

        const double a = d.Dot(d);
        const double b = 2.0 * (w.Dot(d));
        const double c = w.Dot(w) - Raio() * Raio();
        const double delta = b * b - 4.0 * a * c;

        if (delta < 0.0)
        {
            return intersecoes;
        }

        intersecoes.push_back((-b - sqrt(delta)) / 2.0 * a);
        intersecoes.push_back((-b + sqrt(delta)) / 2.0 * a);

        return intersecoes;
    }

private:
    std::string _rotulo;

    TPonto3D _centro;
    double _raio = 0.0;

    TMaterial _material;
};

// ------------------------------------------------------------------------------------------------

class TJanela
{
public:
    TJanela() = default;
    TJanela(
        const TPonto3D& centro,
        double wJanela,
        double hJanela,
        uint16_t wCanvas,
        uint16_t hCanvas
    ) :
        _centro(centro),
        _wJanela(wJanela),
        _hJanela(hJanela),
        _wCanvas(wCanvas),
        _hCanvas(hCanvas)
    {
        _dx = _wJanela / _wCanvas;
        _dy = _hJanela / _hCanvas;
    }

    const TPonto3D& Centro() const { return _centro; }
    double Largura() const { return _wJanela; }
    double Altura() const { return _hJanela; }
    uint16_t LarguraCanvas() const { return _wCanvas; }
    uint16_t AlturaCanvas() const { return _hCanvas; }

    double X(uint16_t coluna) const
    {
        return -0.5 * _wJanela + 0.5 * _dx + coluna * _dx;
    }
    double Y(uint16_t linha) const
    {
        return 0.5 * _hJanela - 0.5 * _dy - linha * _dy;
    }

private:
    TPonto3D _centro;
    double _wJanela = 0.0;
    double _hJanela = 0.0;
    uint16_t _wCanvas = 0;
    uint16_t _hCanvas = 0;

    double _dx = 0.0;
    double _dy = 0.0;
};

// ------------------------------------------------------------------------------------------------

enum class EFormatoImagem { LOG, PPM, BMP };

// ------------------------------------------------------------------------------------------------

std::string Extensao(EFormatoImagem formato)
{
    if (formato == EFormatoImagem::LOG)
    {
        return "log";
    }

    if (formato == EFormatoImagem::PPM)
    {
        return "ppm";
    }

    if (formato == EFormatoImagem::BMP)
    {
        return "bmp";
    }

    return "";
}

// ------------------------------------------------------------------------------------------------

class IArquivoSaida
{
public:
    IArquivoSaida() = default;
    virtual ~IArquivoSaida() = default;

    virtual bool Aberto() const = 0;
    virtual const char* Caminho() const = 0;
    virtual bool Anexa(const TCor& cor) = 0;
    virtual void Flush() = 0;
};

// ------------------------------------------------------------------------------------------------

class TArquivoLOG : public IArquivoSaida
{
public:
    TArquivoLOG() = delete;
    TArquivoLOG(
        const std::string& caminho,
        uint16_t largura,
        uint16_t altura
    )
        : _caminho(caminho),
          _largura(largura),
          _altura(altura)
    {
        _stream.open(_caminho);

        if (_stream.is_open())
        {
            _stream << "W = " << std::to_string(_largura)
                    << "; "
                    << "H = "
                    << std::to_string(altura)
                    << ";\n";
        }
    }
    ~TArquivoLOG()
    {
        _stream.close();
    }

    bool Aberto() const override
    {
        return _stream.is_open();
    }
    const char* Caminho() const override
    {
        return _caminho.c_str();
    }

    bool Anexa(const TCor& cor) override
    {
        const bool aberto = Aberto();

        if (aberto)
        {
            _stream << "X = " << std::to_string(_coordX) << "; "
                    << "Y = " << std::to_string(_coordY) << "; "
                    << "R = " << std::to_string(cor.R()) << "; "
                    << "G = " << std::to_string(cor.G()) << "; "
                    << "B = " << std::to_string(cor.B()) << ";\n";
            
            _coordX++;
            if (_coordX == _largura)
            {
                _coordX = 0;
                _coordY++;
            }
        }
        
        return aberto;
    }

    bool Anexa(const std::string& logMsg)
    {
        const bool aberto = Aberto();

        if (aberto)
        {
            _stream << "    " << logMsg << "\n";
        }
        
        return aberto;
    }

    void Flush() override
    {
        _stream.flush();
    }

private:
    std::string _caminho;
    uint16_t _largura;
    uint16_t _altura;

    std::ofstream _stream;
    uint16_t _coordX = 0;
    uint16_t _coordY = 0;
};

// ------------------------------------------------------------------------------------------------

class TArquivoPPM : public IArquivoSaida
{
public:
    TArquivoPPM() = delete;
    TArquivoPPM(
        const std::string& caminho,
        uint16_t largura,
        uint16_t altura
    )
        : _caminho(caminho),
          _largura(largura),
          _altura(altura)
    {
        _stream.open(_caminho);

        if (_stream.is_open())
        {
            _stream << "P3"
                    << "\n"
                    << std::to_string(_largura)
                    << " "
                    << std::to_string(altura)
                    << "\n"
                    << "255"
                    << "\n";
        }
    }
    ~TArquivoPPM()
    {
        _stream.close();
    }

    bool Aberto() const override
    {
        return _stream.is_open();
    }
    const char* Caminho() const override
    {
        return _caminho.c_str();
    }

    bool Anexa(const TCor& cor) override
    {
        const bool aberto = Aberto();

        if (aberto)
        {
            _stream << std::to_string(cor.R()) << " "
                    << std::to_string(cor.G()) << " "
                    << std::to_string(cor.B()) << "\n";
        }
        
        return aberto;
    }

    void Flush() override
    {
        _stream.flush();
    }

private:
    std::string _caminho;
    uint16_t _largura;
    uint16_t _altura;

    std::ofstream _stream;
};

// ------------------------------------------------------------------------------------------------

class TArquivoBMP : public IArquivoSaida
{
public:
    TArquivoBMP() = delete;
    TArquivoBMP(
        const std::string& caminho,
        uint16_t largura,
        uint16_t altura
    )
        : _caminho(caminho),
          _largura(largura),
          _altura(altura)
    {
        _imagem = new bmp::Bitmap(_largura, _altura);
    }
    ~TArquivoBMP()
    {
        delete _imagem;
    }

    bool Aberto() const override
    {
        return _imagem->width() > 0 && _imagem->height() > 0;
    }
    const char* Caminho() const override
    {
        return _caminho.c_str();
    }

    bool Anexa(const TCor& cor) override
    {
        const bool aberto = Aberto();

        if (aberto)
        {
            _imagem->set(x, y, bmp::Pixel { cor.R(), cor.G(), cor.B() });

            x++;
            if (x == _imagem->width())
            {
                x = 0;
                y++;
            }
        }
        
        return aberto;
    }

    void Flush() override
    {
        _imagem->save(_caminho);
    }

private:
    std::string _caminho;
    uint16_t _largura;
    uint16_t _altura;

    bmp::Bitmap* _imagem = nullptr;
    std::int32_t x = 0;
    std::int32_t y = 0;
};

// ------------------------------------------------------------------------------------------------

namespace FuncoesGerais
{
    uint8_t Trunca(double n)
    {
        return static_cast<uint8_t>(floor(n));
    }

    TCor Vec2Cor(const TVetor3D& v)
    {
        return { Trunca(v.X()), Trunca(v.Y()), Trunca(v.Z()) };
    }

    std::unique_ptr<IArquivoSaida> FabricaArquivo(
        EFormatoImagem formato,
        const std::string& nomeArquivo,
        uint16_t largura,
        uint16_t altura
    )
    {
        IArquivoSaida* arq = nullptr;

        std::string _nomeArquivo = nomeArquivo;
        _nomeArquivo += ".";
        _nomeArquivo += Extensao(formato);

        if (formato == EFormatoImagem::LOG)
        {
            arq = new TArquivoLOG(_nomeArquivo, largura, altura);
        }
        else if (formato == EFormatoImagem::PPM)
        {
            arq = new TArquivoPPM(_nomeArquivo, largura, altura);
        }
        else if (formato == EFormatoImagem::BMP)
        {
            arq = new TArquivoBMP(_nomeArquivo, largura, altura);
        }

        return std::unique_ptr<IArquivoSaida>(arq);
    }
}

// ------------------------------------------------------------------------------------------------

namespace FuncoesGeometricas
{
    TVetor3D Versor(const TPonto3D& pInicio, const TPonto3D& pFim)
    {
        const TVetor3D v = pFim - pInicio;
        const TVetor3D d = v.Normalizado();

        return d;
    }
}

// ------------------------------------------------------------------------------------------------

class TCena3D
{
public:
    TCena3D() = default;
    TCena3D(const TPonto3D& origem, const TJanela& janela) : _janela(janela), _p0(origem) {}

    const TJanela& Janela() const { return _janela; }
    void Janela(const TJanela& janela) { _janela = janela; }

    const TPonto3D& Origem() const { return _p0; }
    void Origem(const TPonto3D& p0) { _p0 = p0; }

    const TCor& BgColor() const { return _bgColor; }
    void BgColor(const TCor& cor) { _bgColor = cor; }

    double IambR() const
    {
        return _iAmbR;
    }
    void IambR(double iAmbR)
    {
        _iAmbR = iAmbR;
    }
    double IambG() const
    {
        return _iAmbG;
    }
    void IambG(double iAmbG)
    {
        _iAmbG = iAmbG;
    }
    double IambB() const
    {
        return _iAmbB;
    }
    void IambB(double iAmbB)
    {
        _iAmbB = iAmbB;
    }

    void Insere(const IEntidade3D& entidade)
    {
        _entidades.push_back(std::unique_ptr<IEntidade3D>(entidade.Copia()));
    }
    void Insere(const IFonteLuminosa& fonte)
    {
        _fontes.push_back(std::unique_ptr<IFonteLuminosa>(fonte.Copia()));
    }

    void Renderizar(IArquivoSaida& arq)
    {
        if (auto arqLog = dynamic_cast<TArquivoLOG*>(&arq))
        {
            _arqLog = arqLog;
        }

        const uint16_t nLinhas = _janela.AlturaCanvas();
        const uint16_t nColunas = _janela.LarguraCanvas();

        const double z = _janela.Centro().Z();
        for (int l = 0; l < nLinhas; l++)
        {
            const double y = _janela.Y(l);

            for (int c = 0; c < nColunas; c++)
            {
                const double x = _janela.X(c);

                arq.Anexa(Cor({ x, y, z }));
            }
        }

        _arqLog = nullptr;
        arq.Flush();
    }

    void Log(const std::string& msg) const
    {
        if (_arqLog != nullptr)
        {
            _arqLog->Anexa(msg);
        }
    }

private:
    TCor Cor(const TPonto3D& p) const
    {
        const TVetor3D d = FuncoesGeometricas::Versor(_p0, p);
        const TRaio3D raio { _p0, d };

        double intersecaoMaisProximaObservador = std::numeric_limits<double>::max();
        IEntidade3D* entidadeMaisProximaObservador = nullptr;

        for (const std::unique_ptr<IEntidade3D>& entidade : _entidades)
        {
            const std::vector<double> intersecoesValidas = IntersecoesValidas(*entidade, raio);
            if (!intersecoesValidas.empty())
            {
                const double intersecaoEntidadeMaisProximaObservador = intersecoesValidas[0];
                if (intersecaoEntidadeMaisProximaObservador < intersecaoMaisProximaObservador)
                {
                    intersecaoMaisProximaObservador = intersecaoEntidadeMaisProximaObservador;
                    entidadeMaisProximaObservador = entidade.get();
                }
            }
        }

        TCor pixel = _bgColor;

        if (entidadeMaisProximaObservador != nullptr)
        {
            if (auto esfera = dynamic_cast<const TEsfera*>(entidadeMaisProximaObservador))
            {
                pixel = Cor(*esfera, raio, intersecaoMaisProximaObservador);
            }
            else if (auto plano = dynamic_cast<const TPlano*>(entidadeMaisProximaObservador))
            {
                pixel = Cor(*plano, raio, intersecaoMaisProximaObservador);
            }
        }

        return pixel;
    }

    TCor Cor(const TEsfera& esfera, const TRaio3D& raio, double ti) const
    {
        const TMaterial& material = esfera.Material();
        const TPonto3D kd { material.KdR(), material.KdG(), material.KdB() };
        const TPonto3D ke { material.KeR(), material.KeG(), material.KeB() };

        const TPonto3D p1 = raio.Ponto(ti);
        const TVetor3D n = TVetor3D { p1 - esfera.Centro() } / esfera.Raio();
        const TVetor3D v = raio.Direcao() * -1.0;

        TVetor3D i;
        for (const std::unique_ptr<IFonteLuminosa>& fonte : _fontes)
        {
            if (auto fontePontual = dynamic_cast<const TFontePontual*>(fonte.get()))
            {
                const TPonto3D pFonte = fontePontual->Posicao();
                const TPonto3D iFonte = fontePontual->Intensidade();
                
                const TVetor3D l = TVetor3D { pFonte - p1 }.Normalizado();
                const TVetor3D r = (n * 2.0 * n.Dot(l)) - l;
                const double fd = std::max(0.0, n.Dot(l));
                const double fe = pow(std::max(v.Dot(r), 0.0), material.M());

                const auto id = TVetor3D(iFonte).Arroba(kd) * fd;
                const auto ie = TVetor3D(iFonte).Arroba(ke) * fe;

                const TVetor3D iCorrente = id + ie;
                i += iCorrente;
            }
        }

        TVetor3D c = i.Arroba(TVetor3D(255.0, 255.0, 255.0));
        c.Clamp(0.0, 255.0);

        return FuncoesGerais::Vec2Cor(c);
    }

    TCor Cor(const TPlano& plano, const TRaio3D& raio, double ti) const
    {
        const TMaterial& material = plano.Material();
        const TPonto3D kd { material.KdR(), material.KdG(), material.KdB() };
        const TPonto3D ke { material.KeR(), material.KeG(), material.KeB() };

        const TPonto3D p1 = raio.Ponto(ti);
        const TVetor3D n = plano.Normal();
        const TVetor3D v = raio.Direcao() * -1.0;

        TVetor3D i;
        for (const std::unique_ptr<IFonteLuminosa>& fonte : _fontes)
        {
            if (auto fontePontual = dynamic_cast<const TFontePontual*>(fonte.get()))
            {
                const TPonto3D pFonte = fontePontual->Posicao();
                const TPonto3D iFonte = fontePontual->Intensidade();
                
                const TVetor3D l = TVetor3D { pFonte - p1 }.Normalizado();
                const TVetor3D r = (n * 2.0 * n.Dot(l)) - l;
                const double fd = std::max(0.0, n.Dot(l));
                const double fe = pow(std::max(v.Dot(r), 0.0), material.M());

                const auto id = TVetor3D(iFonte).Arroba(kd) * fd;
                const auto ie = TVetor3D(iFonte).Arroba(ke) * fe;

                const TVetor3D iCorrente = id + ie;
                i += iCorrente;
            }
        }

        TVetor3D c = i.Arroba(TVetor3D(255.0, 255.0, 255.0));
        c.Clamp(0.0, 255.0);

        return FuncoesGerais::Vec2Cor(c);
    }

    std::vector<double> IntersecoesValidas(
        const IEntidade3D& entidade,
        const TRaio3D& raio
    ) const
    {
        std::vector<double> intersecoesValidas;

        auto TemIntersecao = [&intersecoesValidas](double intersecao)
        {
            for (double intersecaoValida : intersecoesValidas)
            {
                if (intersecao == intersecaoValida)
                {
                    return true;
                }
            }

            return false;
        };

        std::vector<double> intersecoes = entidade.Intersecoes(raio);
        for (double intersecao : intersecoes)
        {
            if (intersecao >= 0.0 && !TemIntersecao(intersecao))
            {
                intersecoesValidas.push_back(intersecao);
            }
        }

        return intersecoesValidas;
    }

    TArquivoLOG* _arqLog = nullptr;

    TJanela _janela;
    TPonto3D _p0; // olho do pintor (origem)

    TCor _bgColor;
    double _iAmbR = 0.0;
    double _iAmbG = 0.0;
    double _iAmbB = 0.0;

    std::vector<std::unique_ptr<IEntidade3D>> _entidades;
    std::vector<std::unique_ptr<IFonteLuminosa>> _fontes;
};

// ------------------------------------------------------------------------------------------------

TEsfera FabricaEsfera()
{
    TMaterial material;
    material.KdR(0.7);
    material.KdG(0.2);
    material.KdB(0.2);
    material.KeR(0.7);
    material.KeG(0.2);
    material.KeB(0.2);
    material.KaR(0.7);
    material.KaG(0.2);
    material.KaB(0.2);
    material.M(10.0);

    TEsfera esfera { { 0.0, 0.0, -100.0 }, 40.0 };
    esfera.Rotulo("ESFERA_1");
    esfera.Material(material);

    return esfera;
}

// ------------------------------------------------------------------------------------------------

TPlano FabricaPlanoChao(const TEsfera& esfera)
{
    TMaterial material;
    material.KdR(0.2);
    material.KdG(0.7);
    material.KdB(0.2);
    material.KeR(0.0);
    material.KeG(0.0);
    material.KeB(0.0);
    material.KaR(0.2);
    material.KaG(0.7);
    material.KaB(0.2);
    material.M(1.0);

    TPlano planoChao { { 0.0, -esfera.Raio(), 0.0 }, { 0.0, 1.0, 0.0 } };
    planoChao.Rotulo("PLANO_CHAO");
    planoChao.Material(material);

    return planoChao;
}

// ------------------------------------------------------------------------------------------------

TPlano FabricaPlanoFundo()
{
    TMaterial material;
    material.KdR(0.3);
    material.KdG(0.3);
    material.KdB(0.7);
    material.KeR(0.0);
    material.KeG(0.0);
    material.KeB(0.0);
    material.KaR(0.3);
    material.KaG(0.3);
    material.KaB(0.7);
    material.M(1.0);

    TPlano planoFundo { { 0.0, 0.0, -200.0 }, { 0.0, 0.0, 1.0 } };
    planoFundo.Rotulo("PLANO_FUNDO");
    planoFundo.Material(material);

    return planoFundo;
}

// ------------------------------------------------------------------------------------------------

TCena3D FabricaCena()
{
    const TPonto3D p0 { 0.0, 0.0, 0.0 };

    const double wJanela = 60.0;
    const double hJanela = 60.0;
    const double dJanela = 30.0;
    const uint16_t wCanvas = 500u;
    const uint16_t hCanvas = 500u;
    const TJanela janela { { 0.0, 0.0, -dJanela }, wJanela, hJanela, wCanvas, hCanvas };

    TCena3D cena { p0, janela };
    cena.BgColor({ 100u, 100u, 100u });
    cena.IambR(0.3);
    cena.IambG(0.3);
    cena.IambB(0.3);

    const TFontePontual fontePontual { { 0.0, 60.0, -30.0 }, { 0.7, 0.7, 0.7 } };
    const TEsfera esfera = FabricaEsfera();
    const TPlano planoChao = FabricaPlanoChao(esfera);
    const TPlano planoFundo = FabricaPlanoFundo();

    cena.Insere(fontePontual);
    cena.Insere(esfera);
    cena.Insere(planoChao);
    cena.Insere(planoFundo);

    return cena;
}

// ------------------------------------------------------------------------------------------------

std::unique_ptr<IArquivoSaida> FabricaArquivo(
    const std::string& nome,
    EFormatoImagem formato,
    const TCena3D& cena
)
{
    const TJanela janela = cena.Janela();

    return FuncoesGerais::FabricaArquivo(
        formato, nome, janela.LarguraCanvas(), janela.AlturaCanvas()
    );
}

// ------------------------------------------------------------------------------------------------

int main()
{
    TCena3D cena = FabricaCena();
    const std::unique_ptr<IArquivoSaida> arq = FabricaArquivo("a", EFormatoImagem::BMP, cena);

    const bool erro = !arq->Aberto();
    if (!erro)
    {
        cena.Renderizar(*arq);
    }

    return erro;
}

// ------------------------------------------------------------------------------------------------
