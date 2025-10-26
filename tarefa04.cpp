// ------------------------------------------------------------------------------------------------
// COMPUTACAO GRAFICA I - TAREFA 04
// ------------------------------------------------------------------------------------------------
// Cena da  Tarefa 03 ampliada
// ------------------------------------------------------------------------------------------------
// Na cena da Tarefa 03, inclua os seguintes objetos
// Cilindro:
// >> Centro da base localizado no centro da esfera
// >> Raio da base igual a im terço do Raio da esfera 
// >> Altura do cilindro igual tres vezes o Raio da esfera
// >> Vetor-direcao do cilindro d_cil = (-1/sqrt(3), 1/sqrt(3), -1/sqrt(3))
// >> Kd = Ke = Ka = ( 0.2, 0.3, 0.8)
// Cone:
// >> Centro da base localizado no centro do topo do cilindro
// >> Raio da base igual a 1.5*Raio da esfera
// >> Altura do cone igual a (1/3)Raio da base
// >> Vetor-direcao do cilindro d_cil = (-1/sqrt(3), 1/sqrt(3), -1/sqrt(3))
// >>  Kd = Ke = Ka = ( 0.8, 0.3, 0.2)
// ------------------------------------------------------------------------------------------------

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <memory>

#define UTILIZA_BITMAP_PLUS_PLUS

#ifdef UTILIZA_BITMAP_PLUS_PLUS
    #include "BitmapPlusPlus.hpp"
#else
    // Fiz a exportacao em bitmap so pra ajudar nas minhas depuracoes,
    // nao vou mandar junto. Se nao tem a biblioteca, usa essa classe
    // dummy que implementa a interface necessaria com stubs, so pro
    // compilador nao reclamar, dai TArquivoBMP nao faz nada nesse caso
    namespace bmp
    {
        struct Pixel { std::uint8_t _1; std::uint8_t _2; std::uint8_t _3; };
        struct Bitmap
        {
            Bitmap(std::int32_t, std::int32_t) {}
            std::int32_t width() const { return 0; };
            std::int32_t height() const { return 0; };
            void set(const std::int32_t, const std::int32_t, const Pixel) {};
            void save(const std::string&) {}
        };
    }
#endif

// ------------------------------------------------------------------------------------------------
#include "agc.hpp"
using TMatriz = agc::matrix<double>;
#define Transposta transpose
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

    TMatriz Vec2Mtx(const TVetor3D& v)
    {
        TMatriz m { 3, 1 };

        m[1][1] = v.X();
        m[2][1] = v.Y();
        m[3][1] = v.Z();

        return m;
    }

    TVetor3D Mtx2Vec(const TMatriz& m)
    {
        return TVetor3D { m[1][1], m[2][1], m[3][1] };
    }

    TMatriz Identidade(int ordem)
    {
        TMatriz I { ordem, ordem };

        for (int i = 1; i <= ordem; i++)
        {
            for (int j = 1; j <= ordem; j++)
            {
                I[i][j] = i == j ? 1.0 : 0.0;
            }
        }

        return I;
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
    virtual TMaterial Material() const = 0;

    virtual TVetor3D Normal(const TPonto3D& p) const = 0;
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

    TMaterial Material() const override
    {
        return _material;
    }
    void Material(const TMaterial& material)
    {
        _material = material;
    }

    const TPonto3D& PontoReferencia() const
    {
        return _p;
    }
    const TVetor3D& Normal() const
    {
        return _n;
    }

    TVetor3D Normal(const TPonto3D& p) const override
    {
        return Normal();
    }

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

protected:
    std::string _rotulo;
    TMaterial _material;

    TPonto3D _p;
    TVetor3D _n;
};

// ------------------------------------------------------------------------------------------------

class TSuperficieCircular : public TPlano
{
public:
    TSuperficieCircular() = delete;
    TSuperficieCircular(const TPonto3D& posicao, const TVetor3D& direcao, double raio)
        : TPlano { posicao, direcao }, _r(raio) {}

    IEntidade3D* Copia() const override
    {
        return new TSuperficieCircular(*this);
    }

    const TPonto3D& PontoReferencia() const
    {
        return _p;
    }

    std::vector<double> Intersecoes(const TRaio3D& raio) const override
    {
        std::vector<double> intersecoes;

        std::vector<double> intersecoesPlano = TPlano::Intersecoes(raio);
        for (double intersecaoPlano : intersecoesPlano)
        {
            if (TVetor3D { raio.Ponto(intersecaoPlano) - _p }.Norma() <= _r)
            {
                intersecoes.push_back(intersecaoPlano);
            }
        }

        return intersecoes;
    }

private:
    double _r;
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

    TMaterial Material() const override
    {
        return _material;
    }
    void Material(const TMaterial& material)
    {
        _material = material;
    }

    const TPonto3D& Centro() const { return _centro; }
    double Raio() const { return _raio; }

    TVetor3D Normal(const TPonto3D& p) const override
    {
        return TVetor3D { p - Centro() } / Raio();
    }

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
    TMaterial _material;

    TPonto3D _centro;
    double _raio = 0.0;
};

// ------------------------------------------------------------------------------------------------

class TSuperficieCilindrica : public IEntidade3D
{
public:
    TSuperficieCilindrica(const TPonto3D& pCentroBase, double raio, double altura, const TVetor3D& direcao)
        : _c(pCentroBase), _r(raio), _h(altura), _d(direcao.Normalizado()), _M(3, 3)
    {
        const TMatriz dc = FuncoesGerais::Vec2Mtx(_d);
        const TMatriz I = FuncoesGerais::Identidade(3);

        _M = I * dc * dc.Transposta();
    }

    IEntidade3D* Copia() const override
    {
        return new TSuperficieCilindrica(*this);
    }

    std::string Rotulo() const override
    {
        return _rotulo;
    }
    void Rotulo(const std::string& rotulo)
    {
        _rotulo = rotulo;
    }

    TMaterial Material() const override
    {
        return _material;
    }
    void Material(const TMaterial& material)
    {
        _material = material;
    }

    TVetor3D Normal(const TPonto3D& p) const override
    {
        return FuncoesGerais::Mtx2Vec(_M * FuncoesGerais::Vec2Mtx(p - _c));
    }
    std::vector<double> Intersecoes(const TRaio3D& raio) const override
    {
        std::vector<double> intersecoes;

        // const TMatriz dr = FuncoesGerais::Vec2Mtx(raio.Direcao());

        // const TMatriz w = FuncoesGerais::Vec2Mtx(raio.Origem() - _c);
        // const TMatriz wT = w.Transposta();
        // // r_ = M__ s_
        // // w_ = p0 - Cb
        // // 0 <= s1_ dot dc_ <= H
        // // nr_ = || ri_ || / Rc
        // // ri_ = M__ si_

        // const double a = TMatriz(dr.Transposta() * _M * dr)[1][1];
        // const double b = TMatriz(2.0 * wT * _M * dr)[1][1];
        // const double c = TMatriz(wT * _M * w)[1][1] - _r * _r;

        const TVetor3D s = raio.Origem() - _c;
        const TVetor3D v = s - _d * s.Dot(_d);
        const TVetor3D w = raio.Direcao() - _d * raio.Direcao().Dot(_d);

        const double a = w.Dot(w);
        const double b = v.Dot(w);
        const double c = v.Dot(v) - _r * _r;

        const double delta = b * b - 4.0 * a * c;

        if (delta < 0.0)
        {
            return intersecoes;
        }

        const double t1 = (-b - sqrt(delta)) / 2.0 * a;
        const TVetor3D s1 = raio.Ponto(t1) - _c;
        const double s1dc = s1.Dot(_d);
        const bool t1Valido = s1dc >= 0.0 && s1dc <= _h;

        const double t2 = (-b + sqrt(delta)) / 2.0 * a;
        const TVetor3D s2 = raio.Ponto(t2) - _c;
        const double s2dc = s2.Dot(_d);
        const bool t2Valido = s2dc >= 0.0 && s2dc <= _h;

        if (t1Valido)
        {
            intersecoes.push_back(t1);
        }
        if (t2Valido)
        {
            intersecoes.push_back(t2);
        }

        return intersecoes;
    }

    const TPonto3D& CentroBase() const
    {
        return _c;
    }
    double Raio() const
    {
        return _r;
    }
    double Altura() const
    {
        return _h;
    }
    const TVetor3D& Direcao() const
    {
        return _d;
    }

private:
    std::string _rotulo;
    TMaterial _material;

    TPonto3D _c;
    double _r;
    double _h;
    TVetor3D _d;

    TVetor3D _n;
    TMatriz _M;
};

// ------------------------------------------------------------------------------------------------

class TCilindro : public IEntidade3D
{
public:
    TCilindro(const TPonto3D& pCentroBase, double raio, double altura, const TVetor3D& direcao)
        : _lateral(pCentroBase, raio, altura, direcao),
          _base(pCentroBase, direcao, raio),
          _topo(pCentroBase + direcao * altura, direcao, raio)
    {
    }

    IEntidade3D* Copia() const override
    {
        return new TCilindro(*this);
    }

    std::string Rotulo() const override
    {
        return _rotulo;
    }
    void Rotulo(const std::string& rotulo)
    {
        _rotulo = rotulo;
    }

    TMaterial Material() const override
    {
        return _material;
    }
    void Material(const TMaterial& material)
    {
        _material = material;
    }

    TVetor3D Normal(const TPonto3D& p) const override
    {
        TVetor3D normal;

        if (_ultima != nullptr)
        {
            normal = _ultima->Normal(p);
        }

        return normal;
    }
    std::vector<double> Intersecoes(const TRaio3D& raio) const override
    {
        std::vector<double> intersecoes;

        const IEntidade3D* entidadeInterceptada = EntidadeInterceptada(raio);
        if (entidadeInterceptada != nullptr)
        {
            intersecoes = entidadeInterceptada->Intersecoes(raio);
        }

        return intersecoes;
    }

    const TPonto3D& CentroBase() const
    {
        return _lateral.CentroBase();
    }
    double Raio() const
    {
        return _lateral.Raio();
    }
    double Altura() const
    {
        return _lateral.Altura();
    }
    const TVetor3D& Direcao() const
    {
        return _lateral.Direcao();
    }

private:
    const IEntidade3D* EntidadeInterceptada(const TRaio3D& raio) const
    {
        const IEntidade3D* entidadeInterceptada = nullptr;

        struct TIntersecoesEntidades
        {
            const IEntidade3D* entidade;
            std::vector<double>* intersecoes;
        };
        std::vector<TIntersecoesEntidades> ies;

        std::vector<double> intersecoesLateral = _lateral.Intersecoes(raio);
        if (!intersecoesLateral.empty())
        {
            ies.push_back({ &_lateral, &intersecoesLateral });
        }

        std::vector<double> intersecoesBase = _base.Intersecoes(raio);
        if (!intersecoesBase.empty())
        {
            ies.push_back({ &_base, &intersecoesBase });
        }

        std::vector<double> intersecoesTopo = _topo.Intersecoes(raio);
        if (!intersecoesTopo.empty())
        {
            ies.push_back({ &_topo, &intersecoesTopo });
        }

        if (!ies.empty())
        {
            std::sort(ies.begin(), ies.end(), [](const auto& i1, const auto& i2)
            {
                return (*i1.intersecoes)[0] < (*i2.intersecoes)[0];
            });

            entidadeInterceptada = ies[0].entidade;
        }

        _ultima = entidadeInterceptada;

        return entidadeInterceptada;
    }

    std::string _rotulo;
    TMaterial _material;

    mutable const IEntidade3D* _ultima = nullptr;

    TSuperficieCilindrica _lateral;
    TSuperficieCircular _base;
    TSuperficieCircular _topo;
};

// ------------------------------------------------------------------------------------------------

class TSuperficieConica : public IEntidade3D
{
public:
    TSuperficieConica(const TPonto3D& pCentroBase, double raioBase, double altura, const TVetor3D& direcao)
        : _c(pCentroBase), _r(raioBase), _h(altura), _d(direcao.Normalizado())
    {
        _v = _c + _d * _h;
    }

    IEntidade3D* Copia() const override
    {
        return new TSuperficieConica(*this);
    }

    std::string Rotulo() const override
    {
        return _rotulo;
    }
    void Rotulo(const std::string& rotulo)
    {
        _rotulo = rotulo;
    }

    TMaterial Material() const override
    {
        return _material;
    }
    void Material(const TMaterial& material)
    {
        _material = material;
    }

    TVetor3D Normal(const TPonto3D& p) const override
    {
        TVetor3D _s = _v - p;
        _s = _s.Normalizado();

        TMatriz s = FuncoesGerais::Vec2Mtx(_s);
        TMatriz M = FuncoesGerais::Identidade(3) - s * s.Transposta();

        return FuncoesGerais::Mtx2Vec(M * FuncoesGerais::Vec2Mtx(_d)).Normalizado();
    }
    std::vector<double> Intersecoes(const TRaio3D& raio) const override
    {
        std::vector<double> intersecoes;

        // const TMatriz dc = FuncoesGerais::Vec2Mtx(_d);
        // const TMatriz I = FuncoesGerais::Identidade(3);
        // const TMatriz M = I * dc * dc.Transposta();
        // const TMatriz dr = FuncoesGerais::Vec2Mtx(raio.Direcao());
        // const TMatriz drT = dr.Transposta();

        // double hr = _h / _r;
        // hr *= hr;

        // const TMatriz Mbarra = dc * dc.Transposta();
        // const TMatriz Masterisco = Mbarra - M * hr;

        // const TMatriz w = FuncoesGerais::Vec2Mtx(raio.Origem() - _c);
        // const TMatriz wT = w.Transposta();

        // const double a = TMatriz(drT * Masterisco * dr)[1][1];
        // const double b = TMatriz(2.0 * wT * Masterisco * dr - 2.0 * _h * drT * dc)[1][1];
        // const double c = TMatriz(wT * Masterisco * w - 2.0 * _h * wT * dc)[1][1] + _h * _h;

        const TPonto3D V = _c + _d * _h;
        const TVetor3D v = V - raio.Origem();
        const double teta = 0.0;
        const double powcosteta = cos(teta) * cos(teta);
        const double a = pow(_d.Dot(raio.Direcao()), 2) - raio.Direcao().Dot(raio.Direcao()) * powcosteta;
        const double b = v.Dot(raio.Direcao()) * powcosteta - v.Dot(_d) * raio.Direcao().Dot(_d);
        const double c = pow(v.Dot(_d), 2) - v.Dot(v) * powcosteta;

        const double delta = b * b - 4.0 * a * c;

        if (delta < 0.0)
        {
            return intersecoes;
        }

        const double t1 = (-b - sqrt(delta)) / 2.0 * a;
        const TVetor3D s1 = raio.Ponto(t1) - _c;
        const double s1dc = s1.Dot(_d);
        const bool t1Valido = s1dc >= 0.0 && s1dc <= _h;

        const double t2 = (-b + sqrt(delta)) / 2.0 * a;
        const TVetor3D s2 = raio.Ponto(t2) - _c;
        const double s2dc = s2.Dot(_d);
        const bool t2Valido = s2dc >= 0.0 && s2dc <= _h;

        if (t1Valido)
        {
            intersecoes.push_back(t1);
        }
        if (t2Valido)
        {
            intersecoes.push_back(t2);
        }

        return intersecoes;
    }

    const TPonto3D& CentroBase() const
    {
        return _c;
    }
    double RaioBase() const
    {
        return _r;
    }
    double Altura() const
    {
        return _h;
    }
    const TVetor3D& Direcao() const
    {
        return _d;
    }

private:
    std::string _rotulo;
    TMaterial _material;

    TPonto3D _c;
    double _r;
    double _h;
    TVetor3D _d;

    TPonto3D _v;
};

// ------------------------------------------------------------------------------------------------

class TCone : public IEntidade3D
{
public:
    TCone(const TPonto3D& pCentroBase, double raio, double altura, const TVetor3D& direcao)
        : _lateral(pCentroBase, raio, altura, direcao),
          _base(pCentroBase, direcao, raio)
    {
    }

    IEntidade3D* Copia() const override
    {
        return new TCone(*this);
    }

    std::string Rotulo() const override
    {
        return _rotulo;
    }
    void Rotulo(const std::string& rotulo)
    {
        _rotulo = rotulo;
    }

    TMaterial Material() const override
    {
        return _material;
    }
    void Material(const TMaterial& material)
    {
        _material = material;
    }

    TVetor3D Normal(const TPonto3D& p) const override
    {
        TVetor3D normal;

        if (_ultima != nullptr)
        {
            normal = _ultima->Normal(p);
        }

        return normal;
    }
    std::vector<double> Intersecoes(const TRaio3D& raio) const override
    {
        std::vector<double> intersecoes;

        const IEntidade3D* entidadeInterceptada = EntidadeInterceptada(raio);
        if (entidadeInterceptada != nullptr)
        {
            intersecoes = entidadeInterceptada->Intersecoes(raio);
        }

        return intersecoes;
    }

    const TPonto3D& CentroBase() const
    {
        return _lateral.CentroBase();
    }
    double RaioBase() const
    {
        return _lateral.RaioBase();
    }
    double Altura() const
    {
        return _lateral.Altura();
    }
    const TVetor3D& Direcao() const
    {
        return _lateral.Direcao();
    }

private:
    const IEntidade3D* EntidadeInterceptada(const TRaio3D& raio) const
    {
        const IEntidade3D* entidadeInterceptada = nullptr;

        struct TIntersecoesEntidades
        {
            const IEntidade3D* entidade;
            std::vector<double>* intersecoes;
        };
        std::vector<TIntersecoesEntidades> ies;

        std::vector<double> intersecoesLateral = _lateral.Intersecoes(raio);
        if (!intersecoesLateral.empty())
        {
            ies.push_back({ &_lateral, &intersecoesLateral });
        }

        std::vector<double> intersecoesBase = _base.Intersecoes(raio);
        if (!intersecoesBase.empty())
        {
            ies.push_back({ &_base, &intersecoesBase });
        }

        if (!ies.empty())
        {
            std::sort(ies.begin(), ies.end(), [](const auto& i1, const auto& i2)
            {
                return (*i1.intersecoes)[0] < (*i2.intersecoes)[0];
            });

            entidadeInterceptada = ies[0].entidade;
        }

        _ultima = entidadeInterceptada;

        return entidadeInterceptada;
    }

    std::string _rotulo;
    TMaterial _material;

    mutable const IEntidade3D* _ultima = nullptr;

    TSuperficieConica _lateral;
    TSuperficieCircular _base;
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

namespace FuncoesGeometricas
{
    TVetor3D Versor(const TPonto3D& pInicio, const TPonto3D& pFim)
    {
        const TVetor3D v = pFim - pInicio;
        const TVetor3D d = v.Normalizado();

        return d;
    }

    std::vector<double> IntersecoesValidas(
        const IEntidade3D& entidade,
        const TRaio3D& raio
    )
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
            const std::vector<double> intersecoesValidas = FuncoesGeometricas::IntersecoesValidas(*entidade, raio);
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
            pixel = Cor(*entidadeMaisProximaObservador, raio, intersecaoMaisProximaObservador);
        }

        return pixel;
    }

    TCor Cor(const IEntidade3D& entidade, const TRaio3D& raio, double ti) const
    {
        const TVetor3D iAmb { IambR(), IambG(), IambB() };

        const TMaterial& material = entidade.Material();
        const TPonto3D kd { material.KdR(), material.KdG(), material.KdB() };
        const TPonto3D ke { material.KeR(), material.KeG(), material.KeB() };
        const TPonto3D ka { material.KaR(), material.KaG(), material.KaB() };

        const TPonto3D pi = raio.Ponto(ti);
        const TVetor3D n = entidade.Normal(pi);
        const TVetor3D v = raio.Direcao() * -1.0;

        const auto ia = iAmb.Arroba(ka);

        TVetor3D i = ia;
        for (const std::unique_ptr<IFonteLuminosa>& fonte : _fontes)
        {
            if (auto fontePontual = dynamic_cast<const TFontePontual*>(fonte.get()))
            {
                const TPonto3D& pFonte = fontePontual->Posicao();
                const TVetor3D& iFonte = fontePontual->Intensidade();
                
                const TVetor3D L = pFonte - pi;
                const TVetor3D l = L.Normalizado();
                const TRaio3D shadowRay { pi, l };

                bool temEntidadeBloqueandoLuz = false;
                for (const std::unique_ptr<IEntidade3D>& outraEntidade : _entidades)
                {
                    if (&entidade != outraEntidade.get())
                    {
                        const std::vector<double> intersecoes = FuncoesGeometricas::IntersecoesValidas(
                            *outraEntidade, shadowRay
                        );
                        
                        if (!intersecoes.empty())
                        {
                            const double intersecaoMaisProxima = intersecoes[0];

                            if (intersecaoMaisProxima < L.Norma())
                            {
                                temEntidadeBloqueandoLuz = true;
                            }
                        }
                    }
                }

                if (!temEntidadeBloqueandoLuz)
                {
                    const TVetor3D r = (n * 2.0 * n.Dot(l)) - l;
                    const double fd = std::max(0.0, n.Dot(l));
                    const double fe = pow(std::max(v.Dot(r), 0.0), material.M());

                    const auto id = iFonte.Arroba(kd) * fd;
                    const auto ie = iFonte.Arroba(ke) * fe;

                    const TVetor3D iCorrente = id + ie;
                    i += iCorrente;
                }
            }
        }

        TVetor3D c = i.Arroba(TVetor3D(255.0, 255.0, 255.0));
        c.Clamp(0.0, 255.0);

        return FuncoesGerais::Vec2Cor(c);
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

TCilindro FabricaCilindro(const TEsfera& ref)
{
    TMaterial material;
    material.KdR(0.2);
    material.KdG(0.3);
    material.KdB(0.8);
    material.KeR(0.2);
    material.KeG(0.3);
    material.KeB(0.8);
    material.KaR(0.2);
    material.KaG(0.3);
    material.KaB(0.8);
    material.M(10.0);

    const TPonto3D& cBaseCilindro = ref.Centro();
    const double rBaseCilindro = ref.Raio() / 3.0;
    const double hCilindro = 3.0 * ref.Raio();
    const double k = 1.0 / sqrt(3.0);
    const TVetor3D dCilindro { -k, k, -k };

    TCilindro cilindro { cBaseCilindro, rBaseCilindro, hCilindro, dCilindro };
    cilindro.Rotulo("CILINDRO_1");
    cilindro.Material(material);

    return cilindro;
}

// ------------------------------------------------------------------------------------------------

TCone FabricaCone(const TEsfera& esferaRef, const TCilindro& cilindroRef)
{
    TMaterial material;
    material.KdR(0.8);
    material.KdG(0.3);
    material.KdB(0.2);
    material.KeR(0.8);
    material.KeG(0.3);
    material.KeB(0.2);
    material.KaR(0.8);
    material.KaG(0.3);
    material.KaB(0.2);
    material.M(10.0);

    const TPonto3D& cBaseCone = cilindroRef.Direcao() * cilindroRef.Altura();
    const double rBaseCone = 1.5 * esferaRef.Raio();
    const double hCone = esferaRef.Raio() / 3.0;
    const TVetor3D& dCone = cilindroRef.Direcao();

    TCone cone { cBaseCone, rBaseCone, hCone, dCone };
    cone.Rotulo("CONE_1");
    cone.Material(material);

    return cone;
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
    const TCilindro cilindro = FabricaCilindro(esfera);
    const TCone cone = FabricaCone(esfera, cilindro);
    const TPlano planoChao = FabricaPlanoChao(esfera);
    const TPlano planoFundo = FabricaPlanoFundo();

    cena.Insere(fontePontual);
    cena.Insere(esfera);
    cena.Insere(cilindro);
    cena.Insere(cone);
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
