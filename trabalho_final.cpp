// ------------------------------------------------------------------------------------------------
// COMPUTACAO GRAFICA I - TRABALHO FINAL
// ------------------------------------------------------------------------------------------------
// 1.     Definicao do cenario 
// 1.1.    Coerencia tematica (Obrigatorio)
// O cenario tem de ter coerencia tematica, isto e, nao pode ser um amontoado de objetos
// aleatoriamente distribuidos.
// 1.2.    Coordenadas do mundo (Obrigatorio) 
// O cenario deve ser montado de forma que todos os objetos estejam no primeiro octante, isto e,
// as coordenadas dos vertices de todos os objetos terao x, y e z positivos.
// 1.3.    Objetos
// 1.3.1. Tipos de objetos (apresentar pelo menos um objeto de cada tipo) (Obrigatorio)
// ·       Esfera
// ·       Cilindro
// ·       Cone
// ·       Malha
// 1.3.2. Materiais (pelo menos quatro materiais distintos) (Obrigatorio)
// 1.3.3. Textura (pelo menos uma textura aplicada) (Obrigatorio)
// 1.4.    Transformacoes
// 1.4.1. Translacao (Obrigatorio)
// 1.4.2. Rotacao
// ·       Em torno de um dos eixos x, y ou z (Obrigatorio)
// ·       Em torno de um eixo arbitrario (Obrigatorio um dos dois metodos)
//              o    Mudanca de sistemas de coordenadas 
//              o    Quaternios
// 1.4.3. Escala (Obrigatorio)
// 1.4.4. Cisalhamento (+ 0.5)
// 1.4.5. Espelho em relacao a um plano arbitrario (+ 0.5)
// 1.5.    Fontes luminosas
// 1.5.1. Pontual (Obrigatorio)
// 1.5.2. Spot (+1.0)
// 1.5.3. Direcional (+0.5)
// 1.5.4. Ambiente (Obrigatorio)
// 2.     Camera
// 2.1.    Permitir a especificacao de (Obrigatorio)
// 2.1.1. Posicao da camera (Eye)
// 2.1.2. Direcionamento de visada (At point)
// 2.1.3. Orientacao da camara em torno do eixo de visada (Up point)
// 2.2.    Parametros adicionais (Obrigatorio)
// 2.2.1. Distancia focal (d)
// 2.2.2. Campo de visao (definir as coordenadas de camera da janela: xmin, xmax, ymin, ymax) 
// 3.     Projecoes
// 3.1.    Perspectiva (Obrigatorio)
// 3.1.1. Alterar os parametros adicionais da camera para  
// 3.1.1.1. aumentar o campo de visao (zoom out) (Obrigatorio)
// 3.1.1.2. diminuir o campo de visao (zoom in)(Obrigatorio)
// 3.1.2. Demonstrar como posicionar a camera para obter
// 3.1.2.1. Perspectiva com um ponto de fuga (+0.5)
// 3.1.2.2. Perspectiva com dois pontos de fuga (+ 0.5)
// 3.1.2.3. Perspectiva com tres ou mais pontos de fuga (+ 0.5)
// 3.2.    Ortografica (+ 0.5)
// 3.3.    Obliqua (+0.5)
// 4.     Sombra (Obrigatorio)
// 5.     Interatividade
// 5.1.    Implementar a funcao de pick (Obrigatorio)
// 5.2.    Uso de interface grafica (Bonus de 0.5 a 1.0)
// 6.     Imagem gerada por ray casting com pelo menos 500 x 500 pixels (Obrigatorio)
// 7.     Bonus de criatividade e beleza (ate 1.0)
// ------------------------------------------------------------------------------------------------

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <memory>
#include <utility>
#include <locale>
#include <codecvt>
#include <iostream>
#include <unordered_map>
#include <thread>
#include <functional>

#ifndef UNICODE
#define UNICODE
#endif
#include <windows.h>
#include <windowsx.h>

#define STB_IMAGE_IMPLEMENTATION
#include "dependencias/stb_image.h"
#include "dependencias/BitmapPlusPlus.hpp"

// ------------------------------------------------------------------------------------------------
#define BUFLEN 256
static bool bp = false;
// ------------------------------------------------------------------------------------------------

namespace ResTbl
{
    static constexpr const char* TEX_MADEIRA  = "recursos/Misc/tex.bmp";
    static constexpr const char* OBJ_GOLDFISH = "recursos/Goldfish/Goldfish_01.obj";
    static constexpr const char* OBJ_SPYRO    = "recursos/Spyro/Spyro.obj";
}

// ------------------------------------------------------------------------------------------------

class Tracer
{
public:
    static void Trace(const std::string& msg)
    {
        if (_isTraceActive)
        {
            if (_traceToFile)
            {
                // Acho que aqui nao precisa setar a flag out explicitamente por ja
                // ser uma ofstream, inclusive funciona sem, mas como nao pesquisei muito
                // e vi que no ctor default da libstdc++ seta, vou deixar assim mesmo
                std::ofstream os(NomeArquivo(), std::ios_base::out | std::ios_base::app);
                if (os.is_open())
                {
                    os << msg << std::endl;
                }
                os.close();
            }
            else
            {
                // Pelo VS, acho que tem que ser assim. Testar depois
                // OutputDebugStringA(msg);
                std::cout << msg << std::endl;
            }
        }
    }

    static void LimpaTrace()
    {
        std::ofstream(NomeArquivo(), std::ios_base::out);
    }

    static void TraceActive(bool traceActive)
    {
        _isTraceActive = traceActive;
    }

    static void TraceToFile(bool traceToFile)
    {
        _traceToFile = traceToFile;
    }

private:
    static constexpr const char* NomeArquivo()
    {
        return "Trace.txt";
    }

    static bool _isTraceActive;
    static bool _traceToFile;
};

bool Tracer::_isTraceActive = true;
bool Tracer::_traceToFile   = true;

// ------------------------------------------------------------------------------------------------

template <typename T, typename S = uint16_t>
class TMatriz
{
public:
    class TLinha // Smart Handle
    {
    public:
        TLinha(T* ptrLinha) : _ptr(ptrLinha) {}

        T& operator[](S coluna)
        {
            return _ptr[coluna - 1];
        }
        const T& operator[](S coluna) const
        {
            return _ptr[coluna - 1];
        }

    private:
        T* _ptr = nullptr;
    };

    TMatriz(S linhas, S colunas)
    {
        _linhas = linhas;
        _colunas = colunas;
        AlocaMem(_linhas, _colunas);
    }

    TMatriz(S linhas, S colunas, const std::vector<std::vector<T>>& elementos)
    {
        _linhas = linhas;
        _colunas = colunas;
        AlocaMem(_linhas, _colunas);

        for (S linha = 0; linha < _linhas; linha++)
        {
            for (S coluna = 0; coluna < _colunas; coluna++)
            {
                _matriz[linha][coluna] = elementos[linha][coluna];
            }
        }
    }

    TMatriz(const TMatriz& outra)
    {
        _linhas = outra._linhas;
        _colunas = outra._colunas;
        AlocaMem(_linhas, _colunas);

        for (S linha = 0; linha < _linhas; linha++)
        {
            for (S coluna = 0; coluna < _colunas; coluna++)
            {
                _matriz[linha][coluna] = outra._matriz[linha][coluna];
            }
        }
    }

    TMatriz(TMatriz&& outra)
    {
        _linhas = outra._linhas;
        outra._linhas = 0;

        _colunas = outra._colunas;
        outra._colunas = 0;

        _matriz = outra._matriz;
        outra._matriz = nullptr;
    }

    TMatriz& operator=(const TMatriz& outra)
    {
        if (&outra != this)
        {
            _linhas = outra._linhas;
            _colunas = outra._colunas;
            AlocaMem(_linhas, _colunas);

            for (S linha = 0; linha < _linhas; linha++)
            {
                for (S coluna = 0; coluna < _colunas; coluna++)
                {
                    _matriz[linha][coluna] = outra._matriz[linha][coluna];
                }
            }
        }

        return *this;
    }

    TMatriz& operator=(TMatriz&& outra)
    {
        if (&outra != this)
        {
            _linhas = outra._linhas;
            outra._linhas = 0;

            _colunas = outra._colunas;
            outra._colunas = 0;

            _matriz = outra._matriz;
            outra._matriz = nullptr;
        }

        return *this;
    }

    ~TMatriz()
    {
        for (S linha = 0; linha < _linhas; linha++)
        {
            delete[] _matriz[linha];
        }

        delete[] _matriz;
    }

    void Inicializa(const T& valorPadrao)
    {
        for (S linha = 0; linha < _linhas; linha++)
        {
            for (S coluna = 0; coluna < _colunas; coluna++)
            {
                _matriz[linha][coluna] = valorPadrao;
            }
        }
    }

    TLinha operator[](S linha)
    {
        return TLinha(_matriz[linha - 1]);
    }
    // caveat: acho que isso permite acesso nao const ao elemento,
    //         porque vai permitir a uma const matrix voltar por copia
    //         um smart handle (nao necessariamente const), que vai
    //         acabar resolvendo no operador nao-const do smart handle,
    //         mas nao quero tratar isso agora, mt coisa ainda pra fazer
    TLinha operator[](S linha) const
    {
        return TLinha(_matriz[linha - 1]);
    }

    S NumeroLinhas() const
    {
        return _linhas;
    }
    S NumeroColunas() const
    {
        return _colunas;
    }

    TMatriz Produto(const TMatriz& outra) const
    {
        if (_colunas != outra._linhas)
        {
            return TMatriz { 0, 0 };
        }

        TMatriz P(_linhas, outra._colunas);
        for (int i = 1; i <= _linhas; i++)
        {
            for (int j = 1; j <= outra._colunas; j++)
            {
                for (int k = 1; k <= _colunas; k++)
                {
                    P[i][j] += (*this)[i][k] * outra[k][j];
                }
            }
        }

        return P;
    }

    bool Inconsistente() const
    {
        return _linhas == 0 || _colunas == 0;
    }

protected:
    void AlocaMem(S linhas, S colunas)
    {
        _matriz = new T*[linhas];
        for (S linha = 0; linha < linhas; linha++)
        {
            // se T for uma classe definida pelo usuario, provavelmente
            // tem um valor default bem definido, mas tipos primitivos
            // seriam inicializados com lixo, entao usamos as chaves para
            // garantir inicializacao semantica (value initialization)
            _matriz[linha] = new T[colunas]{};
        }
    }

    S _linhas;
    S _colunas;

    T** _matriz;
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

    std::string ToString() const
    {
        std::stringstream ss;

        ss << "TCor { "
           << static_cast<int>(_r) << ", "
           << static_cast<int>(_g) << ", "
           << static_cast<int>(_b) << " }";

        return ss.str();
    }

private:
    uint8_t _r = 0;
    uint8_t _g = 0;
    uint8_t _b = 0;
};

// ------------------------------------------------------------------------------------------------

using TImagem = TMatriz<TCor, uint16_t>;

// ------------------------------------------------------------------------------------------------

struct TCoordenadasUV
{
    double u;
    double v;
};

// ------------------------------------------------------------------------------------------------

class TVetor4D
{
public:
    TVetor4D() = default;
    TVetor4D(double x, double y, double z, bool pt = true) : _x(x), _y(y), _z(z), _w(pt) {}

    double X() const { return _x; }
    double Y() const { return _y; }
    double Z() const { return _z; }
    double W() const { return _w; }

    void X(double x) { _x = x; }
    void Y(double y) { _y = y; }
    void Z(double z) { _z = z; }
    void W(double w) { _w = w == 0.0 ? 0 : 1; } // deve ser suficiente pra tratar imprecisoes de PF

    TVetor4D& operator+=(const TVetor4D& outro)
    {
        _x += outro._x;
        _y += outro._y;
        _z += outro._z;
        _w += outro._w;
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
        _w *= k;
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
        _w /= k;
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

    bool EhPonto() const
    {
        return _w == 1;
    }

protected:
    double _x = 0.0;
    double _y = 0.0;
    double _z = 0.0;
    int8_t _w =   1;
};

// ------------------------------------------------------------------------------------------------

class TPonto3D : public TVetor4D
{
public:
    TPonto3D() = default;
    TPonto3D(double x, double y, double z) : TVetor4D(x, y, z) {}
    TPonto3D(const TVetor4D& outro) : TVetor4D(outro.X(), outro.Y(), outro.Z()) {}

    std::string ToString() const
    {
        std::stringstream ss;

        ss << "TPonto3D { " << _x << ", " << _y << ", " << _z << " }";

        return ss.str();
    }
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

    TVetor3D Vetorial(const TVetor3D& outro) const
    {
        return TVetor3D{
            _y * outro._z - _z * outro._y,
            _z * outro._x - _x * outro._z,
            _x * outro._y - _y * outro._x
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

    std::string ToString() const
    {
        std::stringstream ss;

        ss << "TVetor3D { " << _x << ", " << _y << ", " << _z << " }";

        return ss.str();
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

// Eu acho que no mundo Linux o pessoal chama de arquivo tambem, entao daria pra
// ter uma certa liberdade poetica aqui, mas vamos separar os nomes pra ficar mais claro,
// inclusive em relacao as responsabilidades da API de cada um (um dispositivo grafico
// abstrato nao tem um "caminho" no sistema de arquivos, por exemplo)
class IDispositivoSaida
{
public:
    IDispositivoSaida() = default;
    virtual ~IDispositivoSaida() = default;

    virtual bool Anexa(const TCor& cor) = 0;
    virtual void Flush() = 0;
};

// ------------------------------------------------------------------------------------------------

class IArquivoSaida : public IDispositivoSaida
{
public:
    IArquivoSaida() = default;
    virtual ~IArquivoSaida() = default;

    virtual bool Aberto() const = 0;
    virtual const char* Caminho() const = 0;
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

class TWin32Viewport : public IDispositivoSaida
{
public:
    TWin32Viewport() = delete;
    TWin32Viewport(
        HDC hdc,
        std::int32_t x,
        std::int32_t y,
        std::int32_t w,
        std::int32_t h
    ) :
        hdc(hdc), xBeg(x), yBeg(y), w(w), h(h)
    {
    }

    bool Anexa(const TCor& cor) override
    {
        SetPixel(hdc, xBeg + x, yBeg + y, RGB(cor.R(), cor.G(), cor.B()));

        x++;
        if (x == w)
        {
            x = 0;
            y++;
        }

        return true;
    }

    void Flush() override
    {
        // NoOp
    }

private:
    HDC hdc;
    std::int32_t xBeg = 0;
    std::int32_t yBeg = 0;
    std::int32_t w = 0;
    std::int32_t h = 0;

    std::int32_t x = 0;
    std::int32_t y = 0;
};

// ------------------------------------------------------------------------------------------------

class TFrameBuffer : public IDispositivoSaida
{
public:
    TFrameBuffer() = delete;

    TFrameBuffer(uint16_t w, uint16_t h)
    {
        _imagem = new TImagem(h, w);
    }

    TFrameBuffer(const TFrameBuffer& outro)
    {
        _imagem = new TImagem(*outro._imagem);
        _outDevice = outro._outDevice;
    }

    TFrameBuffer(TFrameBuffer&& outro)
    {
        _imagem = outro._imagem;
        outro._imagem = nullptr;

        _outDevice = outro._outDevice;

        _x = outro._x;
        _y = outro._y;
    }

    TFrameBuffer& operator=(const TFrameBuffer& outro)
    {
        if (&outro != this)
        {
            _imagem = new TImagem(*outro._imagem);
            _outDevice = outro._outDevice;

            _x = outro._x;
            _y = outro._y;
        }

        return *this;
    }

    TFrameBuffer& operator=(TFrameBuffer&& outro)
    {
        if (&outro != this)
        {
            _imagem = outro._imagem;
            outro._imagem = nullptr;

            _outDevice = outro._outDevice;

            _x = outro._x;
            _y = outro._y;
        }

        return *this;
    }

    virtual ~TFrameBuffer()
    {
        delete _imagem;
    }

    bool Anexa(const TCor& cor) override
    {
        Pixel(_x, _y, cor);

        _x++;
        if (_x == Largura())
        {
            _x = 0;
            _y++;
        }

        return true;
    }

    // Um frame buffer eh um dispositivo grafico intermediario.
    // Vamos prover uma API para acesso direto, mas, opcionalmente,
    // o usuario pode especificar um dispositivo grafico de saida
    // (como uma TWin32Viewport), de forma que, na chamada de Flush
    // (como em TCena3D::Renderiza), os pixels sejam transferidos
    // de forma gerenciada para esse dispositivo de saida final
    // NOTA: o buffer nao eh "dono" do out device, ou seja, nao eh
    // responsavel por gerenciar sua memoria
    void OutDevice(IDispositivoSaida& outro)
    {
        _outDevice = &outro;
    }

    uint16_t Largura() const
    {
        return _imagem->NumeroColunas();
    }
    uint16_t Altura() const
    {
        return _imagem->NumeroLinhas();
    }

    // O tipo TImagem eh implementado sobre uma classe de matriz subjacente
    // (underlying container) que tem dois detalhes importantes, a saber:
    // 1) a indexacao comeca de 1, nao de 0, como eh comum em programacao,
    //    fiz assim pois, para matrizes, me parece mais natural
    // 2) a ordem de acesso eh mtx[linha][coluna], de forma que temos que
    //    inverter x e y, ja que x seria a coluna, e nao a linha... mas essa
    //    de fato eh a forma canonica de se implementar uma classe de matriz
    const TCor& Pixel(uint16_t x, uint16_t y) const
    {
        return (*_imagem)[y + 1][x + 1];
    }
    void Pixel(uint16_t x, uint16_t y, const TCor& cor)
    {
        (*_imagem)[y + 1][x + 1] = cor;
    }

    void Flush() override
    {
        if (_outDevice != nullptr)
        {
            // Os dispositivos de saida implementam suas logicas de anexar cor
            // de acordo com a ordem de renderizacao da cena, lin > col, isto eh,
            // para cada linha l, para cada coluna c, renderize o pixel (c, l).
            // Portanto, devemos aninhar nosso loop nesse formato para que os pixels
            // sejam escritos nessa ordem que o dispositivo de saida final espera
            for (uint16_t y = 0; y < Altura(); y++)
            {
                for (uint16_t x = 0; x < Largura(); x++)
                {
                    _outDevice->Anexa(Pixel(x, y));
                }
            }

            _outDevice = nullptr;
        }
    }

private:
    TImagem* _imagem = nullptr;
    IDispositivoSaida* _outDevice = nullptr;

    std::uint16_t _x = 0;
    std::uint16_t _y = 0;
};

// ------------------------------------------------------------------------------------------------

namespace FuncoesGerais
{
    bool EstaEm(double val, double min, double max, bool inicioFechado = true, bool fimFechado = true)
    {
        bool maiorQueMin = inicioFechado ? val >= min : val > min;
        bool menorQueMax = fimFechado ? val <= max : val < max;

        return maiorQueMin && menorQueMax;
    }

    uint8_t Trunca(double n)
    {
        return static_cast<uint8_t>(floor(n));
    }

    TCor Vec2Cor(const TVetor3D& v)
    {
        return { Trunca(v.X()), Trunca(v.Y()), Trunca(v.Z()) };
    }
    TVetor3D Cor2Vec(const TCor& c)
    {
        return { c.R() / 255.0, c.G() / 255.0, c.B() / 255.0 };
    }
    TMatriz<double> Vec2Mtx(const TVetor4D& v)
    {
        TMatriz<double> V { 4u, 1u, { { v.X() },
                                      { v.Y() },
                                      { v.Z() },
                                      { v.W() } }};
        return V;
    }
    TVetor4D Mtx2Vec(const TMatriz<double>& mtx)
    {
        TVetor4D v;

        if (mtx.NumeroLinhas() == 4 && mtx.NumeroColunas() == 1)
        {
            v.X(mtx[1][1]);
            v.Y(mtx[2][1]);
            v.Z(mtx[3][1]);
            v.W(mtx[4][1]);
        }

        return v;
    }
    double Deg2Rad(double deg)
    {
        return deg * M_PI / 180.0;
    }
    double Rad2Deg(double rad)
    {
        return 180.0 * rad / M_PI;
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

class TTextura
{
public:
    TTextura()
    {
        _pixels = new TImagem(0, 0);
    }

    TTextura(const TTextura& outra)
    {
        _k = outra._k;
        _pixels = new TImagem(*outra._pixels);
        _caminho = outra._caminho;
    }

    TTextura(TTextura&& outra)
    {
        _k = outra._k;
        _pixels = outra._pixels;
        outra._pixels = nullptr;
        _caminho = outra._caminho;
    }

    TTextura& operator=(const TTextura& outra)
    {
        if (&outra != this)
        {
            _k = outra._k;
            _pixels = new TImagem(*outra._pixels);
            _caminho = outra._caminho;
        }

        return *this;
    }

    TTextura& operator=(TTextura& outra)
    {
        if (&outra != this)
        {
            _k = outra._k;
            _pixels = outra._pixels;
            outra._pixels = nullptr;
            _caminho = outra._caminho;
        }

        return *this;
    }

    virtual ~TTextura()
    {
        delete _pixels;
    }

    void Carrega(const std::string& caminho)
    {
        _caminho = caminho;
        auto ExtensaoEh = [this](const std::string& ext) {
            return _caminho.size() >= ext.size() &&
                   _caminho.compare(_caminho.size() - ext.size(), ext.size(), ext) == 0;
        };

        if (ExtensaoEh(".bmp"))
        {
            CarregaBMP(_caminho);
        }
        else if (ExtensaoEh(".png"))
        {
            CarregaPNG(_caminho);
        }
        else
        {
            _caminho = "";
        }
    }

    TCor Pixel(uint16_t x, uint16_t y) const
    {
        return (*_pixels)[y + 1][x + 1];
    }

    uint16_t Largura() const
    {
        return _pixels->NumeroColunas();
    }
    uint16_t Altura() const
    {
        return _pixels->NumeroLinhas();
    }

    double K() const
    {
        return _k;
    }
    void K(double k)
    {
        _k = k;
    }

    const std::string& Caminho() const
    {
        return _caminho;
    }

private:
    void CarregaBMP(const std::string& caminho)
    {
        bmp::Bitmap img;
        img.load(caminho);

        const uint16_t w = static_cast<uint16_t>(img.width());
        const uint16_t h = static_cast<uint16_t>(img.height());

        delete _pixels;
        _pixels = new TImagem(h, w);

        for (std::size_t x = 0; x < w; x++)
        {
            for (std::size_t y = 0; y < h; y++)
            {
                const bmp::Pixel& pixel = img.get(x, y);
                (*_pixels)[y + 1][x + 1] = { pixel.r, pixel.g, pixel.b };
            }
        }
    }

    void CarregaPNG(const std::string& caminho)
    {
        int w, h, canais;

        unsigned char* data = stbi_load(caminho.c_str(), &w, &h, &canais, 3);

        delete _pixels;
        _pixels = new TImagem(h, w);

        for (int y = 0; y < h; y++)
        {
            for (int x = 0; x < w; x++)
            {
                int idx = (y * w + x) * 3;
                unsigned char r = data[idx];
                unsigned char g = data[idx + 1];
                unsigned char b = data[idx + 2];

                (*_pixels)[y + 1][x + 1] = { r, g, b };
            }
        }

        stbi_image_free(data);
    }

    double _k = 1.0;
    TImagem* _pixels = nullptr;

    std::string _caminho;
};

// ------------------------------------------------------------------------------------------------

class TMaterial
{
public:
    TMaterial() = default;

    TMaterial(const TMaterial& outro)
    {
        CopiaDadosPrimitivos(outro);

        if (outro._tex != nullptr)
        {
            _tex = new TTextura(*outro._tex);
        }
    }

    TMaterial(TMaterial&& outro)
    {
        CopiaDadosPrimitivos(outro);

        _tex = outro._tex;
        outro._tex = nullptr;
    }

    TMaterial& operator=(const TMaterial& outro)
    {
        if (&outro != this)
        {
            CopiaDadosPrimitivos(outro);

            if (outro._tex != nullptr)
            {
                _tex = new TTextura(*outro._tex);
            }
        }

        return *this;
    }

    TMaterial& operator=(TMaterial&& outro)
    {
        if (&outro != this)
        {
            CopiaDadosPrimitivos(outro);

            _tex = outro._tex;
            outro._tex = nullptr;
        }

        return *this;
    }

    virtual ~TMaterial()
    {
        delete _tex;
    }

    double M() const
    {
        return _m;
    }
    void M(double m)
    {
        _m = m;
    }

    double KdR(double u, double v) const
    {
        return _tex != nullptr ? KdTex(u, v).X() : _kdR;
    }
    double KdR() const
    {
        return _kdR;
    }
    void KdR(double r)
    {
        _kdR = r;
    }
    double KdG(double u, double v) const
    {
        return _tex != nullptr ? KdTex(u, v).Y() : _kdG;
    }
    double KdG() const
    {
        return _kdG;
    }
    void KdG(double g)
    {
        _kdG = g;
    }
    double KdB(double u, double v) const
    {
        return _tex != nullptr ? KdTex(u, v).Z() : _kdB;
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

    const TTextura* Textura() const
    {
        return _tex;
    }
    TTextura* Textura()
    {
        return _tex;
    }
    void CarregaTextura(const std::string& caminho)
    {
        _tex = new TTextura;
        _tex->Carrega(caminho);
    }

    std::string ToString() const
    {
        std::string caminhoTextura;
        if (_tex != nullptr)
        {
            caminhoTextura += "\"";
            caminhoTextura += _tex->Caminho();
            caminhoTextura += "\"";
        }
        else
        {
            caminhoTextura = "NULL";
        }

        std::stringstream ss;

        ss << "TMaterial { kd = { "
           << _kdR << ", " << _kdG << ", " << _kdB
           << " }, ke = { "
           << _keR << ", " << _keG << ", " << _keB
           << " }, ka = { "
           << _kaR << ", " << _kaG << ", " << _kaB
           << " }, m = " << _m
           << ", tex = " << caminhoTextura
           << " }";

        return ss.str();
    }

private:
    void CopiaDadosPrimitivos(const TMaterial& outro)
    {
        _m = outro._m;

        _kdR = outro._kdR;
        _kdG = outro._kdG;
        _kdB = outro._kdB;

        _keR = outro._keR;
        _keG = outro._keG;
        _keB = outro._keB;

        _kaR = outro._kaR;
        _kaG = outro._kaG;
        _kaB = outro._kaB;
    }

    TVetor3D KdTex(double u, double v) const
    {
        if (_tex != nullptr)
        {
            auto x = static_cast<uint16_t>(u * (_tex->Largura() - 1));
            auto y = static_cast<uint16_t>((1.0 - v) * (_tex->Altura() - 1));

            return FuncoesGerais::Cor2Vec(_tex->Pixel(x, y));
        }

        return { 0.0, 0.0, 0.0 };
    }

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

    TTextura* _tex = nullptr;
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
    virtual TMaterial Material(const TRaio3D& raio) const = 0;
    virtual TCoordenadasUV CoordenadasUV(const TPonto3D& p, const TRaio3D& raio) const = 0;

    virtual TVetor3D Normal(const TPonto3D& p, const TRaio3D& raio) const = 0;
    virtual std::vector<double> Intersecoes(const TRaio3D& raio) const = 0;

    virtual void Transforma(const TMatriz<double>& matriz4x4) = 0;
};

// ------------------------------------------------------------------------------------------------

namespace TransformacoesLineares
{
    // NOTA: as transformacoes sobre entidades 3D foram colocadas neste namespace
    //       por simplicidade, mas o ideal provavelmente seria ter uma classe de
    //       Transformer (transformador de entidades 3d), que empilhasse as operacoes
    //       a serem aplicadas, otimizando o processamento ao fazer somente 1 produto
    //       de matrizes (pela matriz combinada das transformacoes)

    // A rigor, translacao nao eh uma TL, mas usamos coordenadas homogeneas
    // a fim de poder trata-la como tal
    void Translacao(IEntidade3D& entidade, const TPonto3D& t)
    {
        TMatriz<double> T(4u, 4u, { { 1.0, 0.0, 0.0, t.X() },
                                    { 0.0, 1.0, 0.0, t.Y() },
                                    { 0.0, 0.0, 1.0, t.Z() },
                                    { 0.0, 0.0, 0.0,   1.0 } });
        entidade.Transforma(T);
    }

    void RotacaoEmTornoDeX(IEntidade3D& entidade, double teta)
    {
        TMatriz<double> Rx(4u, 4u, { { 1.0,       0.0,        0.0, 0.0 },
                                     { 0.0, cos(teta), -sin(teta), 0.0 },
                                     { 0.0, sin(teta),  cos(teta), 0.0 },
                                     { 0.0,       0.0,        0.0, 1.0 } });
        entidade.Transforma(Rx);
    }

    void RotacaoEmTornoDeY(IEntidade3D& entidade, double teta)
    {
        TMatriz<double> Ry(4u, 4u, { {  cos(teta), 0.0, sin(teta), 0.0 },
                                     {        0.0, 1.0,       0.0, 0.0 },
                                     { -sin(teta), 0.0, cos(teta), 0.0 },
                                     {        0.0, 0.0,       0.0, 1.0 } });
        entidade.Transforma(Ry);
    }

    void RotacaoEmTornoDeZ(IEntidade3D& entidade, double teta)
    {
        TMatriz<double> Rz(4u, 4u, { { cos(teta), -sin(teta), 0.0, 0.0 },
                                     { sin(teta),  cos(teta), 0.0, 0.0 },
                                     {       0.0,        0.0, 1.0, 0.0 },
                                     {       0.0,        0.0, 0.0, 1.0 } });
        entidade.Transforma(Rz);
    }

    void Rotacao(IEntidade3D& entidade, const TVetor3D& u, double teta)
    {
        // ToDo
    }

    void Escala(IEntidade3D& entidade, const TPonto3D& s)
    {
        TMatriz<double> S(4u, 4u, { { s.X(),   0.0,   0.0, 0.0 },
                                    {   0.0, s.Y(),   0.0, 0.0 },
                                    {   0.0,   0.0, s.Z(), 0.0 },
                                    {   0.0,   0.0,   0.0, 1.0 } });
        entidade.Transforma(S);
    }

    void Cisalhamento(IEntidade3D& entidade, double hXY, double hXZ, double hYX, double hYZ, double hZX, double hZY)
    {
        TMatriz<double> H(4u, 4u, { { 1.0, hXY, hXZ, 0.0 },
                                    { hYX, 1.0, hYZ, 0.0 },
                                    { hZX, hZY, 1.0, 0.0 },
                                    { 0.0, 0.0, 0.0, 1.0 } });
        entidade.Transforma(H);
    }

    void ReflexaoEmTornoDeXY(IEntidade3D& entidade)
    {
        TMatriz<double> R(4u, 4u, { { 1.0, 0.0,  0.0, 0.0 },
                                    { 0.0, 1.0,  0.0, 0.0 },
                                    { 0.0, 0.0, -1.0, 0.0 },
                                    { 0.0, 0.0,  0.0, 1.0 } });
        entidade.Transforma(R);
    }

    void ReflexaoEmTornoDeXZ(IEntidade3D& entidade)
    {
        TMatriz<double> R(4u, 4u, { { 1.0,  0.0, 0.0, 0.0 },
                                    { 0.0, -1.0, 0.0, 0.0 },
                                    { 0.0,  0.0, 1.0, 0.0 },
                                    { 0.0,  0.0, 0.0, 1.0 } });
        entidade.Transforma(R);
    }

    void ReflexaoEmTornoDeYZ(IEntidade3D& entidade)
    {
        TMatriz<double> R(4u, 4u, { { -1.0, 0.0, 0.0, 0.0 },
                                    {  0.0, 1.0, 0.0, 0.0 },
                                    {  0.0, 0.0, 1.0, 0.0 },
                                    {  0.0, 0.0, 0.0, 1.0 } });
        entidade.Transforma(R);
    }

    void Reflexao(IEntidade3D& entidade, const TVetor3D& n, const TPonto3D& Q)
    {
        const TVetor3D _n = n.Normalizado();
        const double a = _n.X();
        const double b = _n.Y();
        const double c = _n.Z();

        TMatriz<double> Tvai(4u, 4u, { { 1.0, 0.0, 0.0, -Q.X() },
                                       { 0.0, 1.0, 0.0, -Q.Y() },
                                       { 0.0, 0.0, 1.0, -Q.Z() },
                                       { 0.0, 0.0, 0.0,    1.0 } });

        TMatriz<double> Rorigem(4u, 4u, { { 1.0 - 2.0 * a * a,      -2.0 * a * b,      -2.0 * a * c, 0.0 },
                                          {      -2.0 * b * a, 1.0 - 2.0 * b * b,      -2.0 * b * c, 0.0 },
                                          {      -2.0 * c * a,      -2.0 * c * b, 1.0 - 2.0 * c * c, 0.0 },
                                          {               0.0,               0.0,               0.0, 1.0 } });

        TMatriz<double> Tvolta(4u, 4u, { { 1.0, 0.0, 0.0, Q.X() },
                                         { 0.0, 1.0, 0.0, Q.Y() },
                                         { 0.0, 0.0, 1.0, Q.Z() },
                                         { 0.0, 0.0, 0.0,   1.0 } });

        TMatriz<double> R = Tvolta.Produto(Rorigem.Produto(Tvai));

        entidade.Transforma(R);
    }

    // Helpers de conveniencia pra nao ter que tratar com TVetor4D diretamente
    bool Transforma(TPonto3D& p, const TMatriz<double>& matriz4x4)
    {
        bool transformou = false;

        const TMatriz<double> P = FuncoesGerais::Vec2Mtx(p);
        const TMatriz<double> P_ = matriz4x4.Produto(P);
        const TVetor4D p_ = FuncoesGerais::Mtx2Vec(P_);

        if (p_.EhPonto())
        {
            p = p_;
            transformou = true;
        }

        return transformou;
    }

    bool Transforma(TVetor3D& v, const TMatriz<double>& matriz4x4)
    {
        bool transformou = false;

        const TMatriz<double> V = FuncoesGerais::Vec2Mtx(v);
        const TMatriz<double> V_ = matriz4x4.Produto(V);
        const TVetor4D v_ = FuncoesGerais::Mtx2Vec(V_);

        if (!v_.EhPonto())
        {
            v = v_;
            transformou = true;
        }

        return transformou;
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

    TVetor3D Normal(const TPonto3D& p1, const TPonto3D& p2, const TPonto3D& p3)
    {
        const TVetor3D r1 = p2 - p1;
        const TVetor3D r2 = p3 - p1;

        return r1.Vetorial(r2);
    }

    double ProdutoMisto(const TVetor3D& v1, const TVetor3D& v2, const TVetor3D& v3)
    {
        return v1.Dot(v2.Vetorial(v3));
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

class TEntidadeComposta : public IEntidade3D
{
public:
    TEntidadeComposta() = default;

    TEntidadeComposta(const TEntidadeComposta& outra)
    {
        for (const std::unique_ptr<IEntidade3D>& entidade : outra._entidades)
        {
            _entidades.push_back(std::unique_ptr<IEntidade3D>(entidade->Copia()));
        }

        _rotulo = outra._rotulo;
    }

    IEntidade3D* Copia() const override
    {
        return new TEntidadeComposta(*this);
    }

    std::string Rotulo() const override
    {
        return _rotulo;
    }
    std::string Rotulo(const TRaio3D& raio) const
    {
        std::string rotulo;

        const IEntidade3D* entidade = EntidadeInterceptada(raio);
        if (entidade != nullptr)
        {
            rotulo = entidade->Rotulo();
        }

        return rotulo;
    }
    void Rotulo(const std::string& rotulo)
    {
        _rotulo = rotulo;
    }

    TMaterial Material(const TRaio3D& raio) const override
    {
        TMaterial material;

        const IEntidade3D* entidade = EntidadeInterceptada(raio);
        if (entidade != nullptr)
        {
            material = entidade->Material(raio);
        }

        return material;
    }

    TCoordenadasUV CoordenadasUV(const TPonto3D& p, const TRaio3D& raio) const override
    {
        TCoordenadasUV uv;

        const IEntidade3D* entidade = EntidadeInterceptada(raio);
        if (entidade != nullptr)
        {
            uv = entidade->CoordenadasUV(p, raio);
        }

        return uv;
    }

    TVetor3D Normal(const TPonto3D& p, const TRaio3D& raio) const override
    {
        TVetor3D normal;

        const IEntidade3D* entidade = EntidadeInterceptada(raio);
        if (entidade != nullptr)
        {
            normal = entidade->Normal(p, raio);
        }

        return normal;
    }
    std::vector<double> Intersecoes(const TRaio3D& raio) const override
    {
        std::vector<double> intersecoes;

        const IEntidade3D* entidade = EntidadeInterceptada(raio);
        if (entidade != nullptr)
        {
            intersecoes = entidade->Intersecoes(raio);
        }

        return intersecoes;
    }

    void Insere(const IEntidade3D& entidade)
    {
        _entidades.push_back(std::unique_ptr<IEntidade3D>(entidade.Copia()));
    }

    void Transforma(const TMatriz<double>& matriz4x4) override
    {
        for (std::unique_ptr<IEntidade3D>& entidade : _entidades)
        {
            entidade->Transforma(matriz4x4);
        }
    }

protected:
    IEntidade3D* EntidadeInterceptada(const TRaio3D& raio) const
    {
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

        return entidadeMaisProximaObservador;
    }

    std::string _rotulo;

    std::vector<std::unique_ptr<IEntidade3D>> _entidades;
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

    TMaterial Material(const TRaio3D&) const override
    {
        return _material;
    }
    void Material(const TMaterial& material)
    {
        _material = material;
    }

    TCoordenadasUV CoordenadasUV(const TPonto3D& p, const TRaio3D& raio) const override
    {
        double u = 0.0;
        double v = 0.0;

        if (const TTextura* tex = _material.Textura())
        {
            u = p.X() / tex->K();
            u = u - floor(u);

            v = p.Z() / tex->K();
            v = v - floor(v);
        }

        return { u, v };
    }

    const TPonto3D& PontoReferencia() const
    {
        return _p;
    }
    const TVetor3D& Normal() const
    {
        return _n;
    }

    TVetor3D Normal(const TPonto3D&, const TRaio3D&) const override
    {
        return Normal();
    }

    std::vector<double> Intersecoes(const TRaio3D& raio) const override
    {
        std::vector<double> intersecoes;

        const double dn = raio.Direcao().Dot(_n);
        if (dn != 0.0)
        {
            const double t = -(_n.Dot(raio.Origem() - _p) / dn);
            intersecoes.push_back(t);
        }

        return intersecoes;
    }

    void Transforma(const TMatriz<double>& matriz4x4) override
    {
        // ToDo
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

    void Transforma(const TMatriz<double>& matriz4x4) override
    {
        // ToDo
    }

private:
    double _r;
};

// ------------------------------------------------------------------------------------------------

class TSuperficieTriangular : public TPlano
{
public:
    TSuperficieTriangular() = delete;
    TSuperficieTriangular(const TPonto3D& p1, const TPonto3D& p2, const TPonto3D& p3)
        : _p1(p1), _p2(p2), _p3(p3), TPlano { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } }
    {
        AtualizaDados();
    }

    IEntidade3D* Copia() const override
    {
        return new TSuperficieTriangular(*this);
    }

    TCoordenadasUV CoordenadasUV(const TPonto3D& p, const TRaio3D& raio) const override
    {
        const double normaN = _N.Norma();
        const double areaABC = FuncoesGeometricas::ProdutoMisto(_p2-_p1, _p3-_p1, _N) / normaN;
        const double volumeABC = normaN * areaABC;

        const double c1 = FuncoesGeometricas::ProdutoMisto(_p2 - p, _p3 - p, _N) / volumeABC;
        const double c2 = FuncoesGeometricas::ProdutoMisto(_p3 - p, _p1 - p, _N) / volumeABC;
        const double c3 = 1.0 - c1 - c2;

        const double u = c1 * _uv1.u + c2 * _uv2.u + c3 * _uv3.u;
        const double v = c1 * _uv1.v + c2 * _uv2.v + c3 * _uv3.v;

        return { u, v };
    }

    const TPonto3D& P1() const
    {
        return _p1;
    }
    const TPonto3D& P2() const
    {
        return _p2;
    }
    const TPonto3D& P3() const
    {
        return _p3;
    }

    void CoordenadasUV1(TCoordenadasUV uv)
    {
        _uv1 = uv;
    }
    void CoordenadasUV2(TCoordenadasUV uv)
    {
        _uv2 = uv;
    }
    void CoordenadasUV3(TCoordenadasUV uv)
    {
        _uv3 = uv;
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
            const TPonto3D pI = raio.Ponto(intersecaoPlano);

            const TVetor3D s1 = _p1 - pI;
            const TVetor3D s2 = _p2 - pI;
            const TVetor3D s3 = _p3 - pI;

            const double normaN = _N.Norma();
            const double c1 = FuncoesGeometricas::ProdutoMisto(_n, s3, s1) / normaN;
            const double c2 = FuncoesGeometricas::ProdutoMisto(_n, s1, s2) / normaN;
            const double c3 = 1.0 - c1 - c2;
            const bool intersecaoDentro = c1 > 0.0 && c2 > 0.0 && c3 > 0.0;

            if (intersecaoDentro)
            {
                intersecoes.push_back(intersecaoPlano);
            }
        }

        return intersecoes;
    }

    void Transforma(const TMatriz<double>& matriz4x4) override
    {
        // ToDo
    }

private:
    void AtualizaDados()
    {
        _p = _p1;
        _N = FuncoesGeometricas::Normal(_p1, _p2, _p3);
        _n = _N.Normalizado();
    }

    TPonto3D _p1;
    TPonto3D _p2;
    TPonto3D _p3;

    TCoordenadasUV _uv1;
    TCoordenadasUV _uv2;
    TCoordenadasUV _uv3;

    TVetor3D _N;
};

// ------------------------------------------------------------------------------------------------

class TSuperficieRetangular : public TEntidadeComposta
{
public:
    TSuperficieRetangular(
        const TPonto3D& p1,
        const TPonto3D& p2,
        const TPonto3D& p3,
        const TPonto3D& p4
    )
    {
        Insere(TSuperficieTriangular { p1, p2, p3 });
        Insere(TSuperficieTriangular { p1, p3, p4 });

        AtualizaReferencias();
    }

    TSuperficieRetangular(const TSuperficieRetangular& outra) : TEntidadeComposta(outra)
    {
        AtualizaReferencias();
    }

    IEntidade3D* Copia() const override
    {
        return new TSuperficieRetangular(*this);
    }

    void Material(const TMaterial& material)
    {
        _t1->Material(material);
        _t2->Material(material);
    }

private:
    void AtualizaReferencias()
    {
        _t1 = static_cast<TSuperficieTriangular*>(_entidades[0].get());
        _t2 = static_cast<TSuperficieTriangular*>(_entidades[1].get());
    }

    TSuperficieTriangular* _t1 = nullptr;
    TSuperficieTriangular* _t2 = nullptr;
};

// ------------------------------------------------------------------------------------------------

class TCubo : public TEntidadeComposta
{
public:
    TCubo(const TPonto3D& pCentro, double aresta) : _a(aresta), _sa(0.5 * aresta)
    {
        const TPonto3D pFaceInf1 { pCentro.X() - _sa, pCentro.Y(), pCentro.Z() + _sa };
        const TPonto3D pFaceInf2 { pCentro.X() + _sa, pCentro.Y(), pCentro.Z() + _sa };
        const TPonto3D pFaceInf3 { pCentro.X() + _sa, pCentro.Y(), pCentro.Z() - _sa };
        const TPonto3D pFaceInf4 { pCentro.X() - _sa, pCentro.Y(), pCentro.Z() - _sa };
        
        const TPonto3D pFaceSup1 { pFaceInf1.X(), pFaceInf1.Y() + _a, pFaceInf1.Z() };
        const TPonto3D pFaceSup2 { pFaceInf2.X(), pFaceInf2.Y() + _a, pFaceInf2.Z() };
        const TPonto3D pFaceSup3 { pFaceInf3.X(), pFaceInf3.Y() + _a, pFaceInf3.Z() };
        const TPonto3D pFaceSup4 { pFaceInf4.X(), pFaceInf4.Y() + _a, pFaceInf4.Z() };

        const TSuperficieRetangular faceInferior { pFaceInf1, pFaceInf2, pFaceInf3, pFaceInf4 };
        const TSuperficieRetangular faceSuperior { pFaceSup1, pFaceSup2, pFaceSup3, pFaceSup4 };
        const TSuperficieRetangular faceLateralEsquerda { pFaceInf1, pFaceSup1, pFaceSup4, pFaceInf4 };
        const TSuperficieRetangular faceLateralDireita { pFaceInf2, pFaceInf3, pFaceSup3, pFaceSup2 };
        const TSuperficieRetangular faceTraseira { pFaceInf3, pFaceInf4, pFaceSup4, pFaceSup3 };
        const TSuperficieRetangular faceFrontal { pFaceInf1, pFaceInf2, pFaceSup2, pFaceSup1 };

        Insere(faceInferior);
        Insere(faceSuperior);
        Insere(faceLateralEsquerda);
        Insere(faceLateralDireita);
        Insere(faceTraseira);
        Insere(faceFrontal);

        AtualizaReferencias();
    }

    TCubo(const TCubo& outro) : TEntidadeComposta(outro)
    {
        AtualizaReferencias();
    }

    IEntidade3D* Copia() const override
    {
        return new TCubo(*this);
    }

    void Material(const TMaterial& material)
    {
        _fInferior->Material(material);
        _fSuperior->Material(material);
        _fLateralEsq->Material(material);
        _fLateralDir->Material(material);
        _fTraseira->Material(material);
        _fFrontal->Material(material);
    }

private:
    void AtualizaReferencias()
    {
        _fInferior = static_cast<TSuperficieRetangular*>(_entidades[0].get());
        _fSuperior = static_cast<TSuperficieRetangular*>(_entidades[1].get());
        _fLateralEsq = static_cast<TSuperficieRetangular*>(_entidades[2].get());
        _fLateralDir = static_cast<TSuperficieRetangular*>(_entidades[3].get());
        _fTraseira = static_cast<TSuperficieRetangular*>(_entidades[4].get());
        _fFrontal = static_cast<TSuperficieRetangular*>(_entidades[5].get());
    }

    double _a = 0.0;
    double _sa = 0.0;

    TSuperficieRetangular* _fInferior = nullptr;
    TSuperficieRetangular* _fSuperior = nullptr;
    TSuperficieRetangular* _fLateralEsq = nullptr;
    TSuperficieRetangular* _fLateralDir = nullptr;
    TSuperficieRetangular* _fTraseira = nullptr;
    TSuperficieRetangular* _fFrontal = nullptr;
};

// ------------------------------------------------------------------------------------------------

class TMalha3D : public TEntidadeComposta
{
public:
    TMalha3D() = default;

    TMalha3D(const TMalha3D& outra) : TEntidadeComposta(outra)
    {
        _deveRotularTriangulosInseridosAutomaticamente = outra._deveRotularTriangulosInseridosAutomaticamente;
    }

    IEntidade3D* Copia() const override
    {
        return new TMalha3D(*this);
    }

    void Visita(std::function<void(const TSuperficieTriangular&)> FuncaoVisita) const
    {
        for (const std::unique_ptr<IEntidade3D>& entidade : _entidades)
        {
            // possivel quirk de sptrs? o metodo eh const, o for ranged base exige
            // que o unique_ptr seja const por isso, mas isso nao garante que o obj
            // para o qual ele aponta eh const, diferente da semantica de raw ptrs,
            // onde um const T* aponta para um objeto t const. Mas aqui, usamos
            // o smart pointer const, e acessamos tranquilamente t como nao-const,
            // ou seja, const unique_ptr<T> nao faz o mesmo efeito de const T*
            // estudar melhor depois
            auto& triangulo = static_cast<TSuperficieTriangular&>(*entidade.get());
            FuncaoVisita(triangulo);
        }
    }

    void Transforma(std::function<void(TSuperficieTriangular&)> FuncaoTransformacao)
    {
        for (std::unique_ptr<IEntidade3D>& entidade : _entidades)
        {
            auto& triangulo = static_cast<TSuperficieTriangular&>(*entidade.get());
            FuncaoTransformacao(triangulo);
        }
    }

    void Adiciona(const TSuperficieTriangular& triangulo)
    {
        TEntidadeComposta::Insere(triangulo);

        if (_deveRotularTriangulosInseridosAutomaticamente)
        {
            auto& trianguloInserido = static_cast<TSuperficieTriangular&>(*_entidades.back().get());
            if (Rotulo() != "" && trianguloInserido.Rotulo() == "")
            {
                std::stringstream ss;
                ss << Rotulo() << "_TRIANGULO_" << _entidades.size();
                trianguloInserido.Rotulo(ss.str());
            }
        }
    }

    TPonto3D Centro() const
    {
        if (_entidades.empty())
        {
            return TPonto3D { 0.0, 0.0, 0.0 };
        }

        auto& trianguloRef = static_cast<TSuperficieTriangular&>(*_entidades.front().get());
        TPonto3D bbMin = trianguloRef.P1();
        TPonto3D bbMax = trianguloRef.P1();

        auto ComputaVertice = [&bbMin, &bbMax](const TPonto3D& v)
        {
            bbMin = TPonto3D(
                std::min(bbMin.X(), v.X()),
                std::min(bbMin.Y(), v.Y()),
                std::min(bbMin.Z(), v.Z())
            );

            bbMax = TPonto3D(
                std::max(bbMax.X(), v.X()),
                std::max(bbMax.Y(), v.Y()),
                std::max(bbMax.Z(), v.Z())
            );
        };

        Visita([ComputaVertice](const TSuperficieTriangular& t) {
            ComputaVertice(t.P1());
            ComputaVertice(t.P2());
            ComputaVertice(t.P3());
        });

        return TPonto3D {
            (bbMin.X() + bbMax.X()) * 0.5,
            (bbMin.Y() + bbMax.Y()) * 0.5,
            (bbMin.Z() + bbMax.Z()) * 0.5
        };
    }

    bool DeveRotularTriangulosInseridosAutomaticamente() const
    {
        return _deveRotularTriangulosInseridosAutomaticamente;
    }
    void DeveRotularTriangulosInseridosAutomaticamente(bool deveRotularAutomaticamente)
    {
        _deveRotularTriangulosInseridosAutomaticamente = deveRotularAutomaticamente;
    }

private:
    bool _deveRotularTriangulosInseridosAutomaticamente = false;
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

    TMaterial Material(const TRaio3D&) const override
    {
        return _material;
    }
    void Material(const TMaterial& material)
    {
        _material = material;
    }

    TCoordenadasUV CoordenadasUV(const TPonto3D& p, const TRaio3D& raio) const override
    {
        return { 0.0, 0.0 };
    }

    const TPonto3D& Centro() const { return _centro; }
    double Raio() const { return _raio; }

    TVetor3D Normal(const TPonto3D& p, const TRaio3D&) const override
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

    void Transforma(const TMatriz<double>& matriz4x4) override
    {
        // ToDo
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
        : _c(pCentroBase), _r(raio), _h(altura), _d(direcao.Normalizado()) {}

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

    TMaterial Material(const TRaio3D&) const override
    {
        return _material;
    }
    void Material(const TMaterial& material)
    {
        _material = material;
    }

    TCoordenadasUV CoordenadasUV(const TPonto3D& p, const TRaio3D& raio) const override
    {
        return { 0.0, 0.0 };
    }

    // nos dois metodos abaixo, as notas de aula estavam bem diferentes
    // das fotos da lousa. inicialmente tentei implementar pelas fotos,
    // mas ocorreram diversos problemas, entao fui pelas notas, que
    // pareceu funcionar com alguns ajustes. depois avaliar o que errei
    // ao fazer a versao da lousa (p.s: aparentemente era bem mais
    // pesada, pelo uso ostensivo de matrizes)
    TVetor3D Normal(const TPonto3D& p, const TRaio3D&) const override
    {
        const TVetor3D s = p - _c;
        const TVetor3D s_ = _d * s.Dot(_d);

        return TVetor3D { s - s_ }.Normalizado();
    }
    std::vector<double> Intersecoes(const TRaio3D& raio) const override
    {
        std::vector<double> intersecoes;

        const TVetor3D& dc = _d;
        const TVetor3D& dr = raio.Direcao();

        const TVetor3D aux = raio.Origem() - _c;
        const TVetor3D v = aux - dc * aux.Dot(dc);
        const TVetor3D w = dr - dc * dr.Dot(dc);

        const double a = w.Dot(w);
        const double b = 2.0 * v.Dot(w);
        const double c = v.Dot(v) - _r * _r;

        // talvez fosse bom avaliar "a" como fabs(a) < 1e-12
        const double delta = b * b - 4.0 * a * c;
        if (a == 0.0 || delta < 0.0)
        {
            return intersecoes;
        }

        // passar para o solver generico depois
        const double t1 = (-b - sqrt(delta)) / (2.0 * a);
        const double t2 = (-b + sqrt(delta)) / (2.0 * a);
        const TPonto3D p1 = raio.Ponto(t1);
        const TPonto3D p2 = raio.Ponto(t2);
        const TVetor3D s1 = p1 - _c;
        const TVetor3D s2 = p2 - _c;
        const double alturaT1 = s1.Dot(dc);
        const double alturaT2 = s2.Dot(dc);
        const bool t1Valido = FuncoesGerais::EstaEm(alturaT1, 0.0, _h);
        const bool t2Valido = FuncoesGerais::EstaEm(alturaT2, 0.0, _h);

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

    void Transforma(const TMatriz<double>& matriz4x4) override
    {
        // ToDo
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
};

// ------------------------------------------------------------------------------------------------

class TCilindro : public TEntidadeComposta
{
public:
    TCilindro(const TPonto3D& pCentroBase, double raio, double altura, const TVetor3D& direcao)
    {
        Insere(TSuperficieCilindrica { pCentroBase, raio, altura, direcao });
        Insere(TSuperficieCircular { pCentroBase, direcao, raio });
        Insere(TSuperficieCircular { pCentroBase + direcao * altura, direcao, raio });

        _lateral = static_cast<TSuperficieCilindrica*>(_entidades[0].get());
        _base = static_cast<TSuperficieCircular*>(_entidades[1].get());
        _topo = static_cast<TSuperficieCircular*>(_entidades[2].get());
    }

    void Material(const TMaterial& material)
    {
        _lateral->Material(material);
        _base->Material(material);
        _topo->Material(material);
    }

    const TPonto3D& CentroBase() const
    {
        return _lateral->CentroBase();
    }
    double Raio() const
    {
        return _lateral->Raio();
    }
    double Altura() const
    {
        return _lateral->Altura();
    }
    const TVetor3D& Direcao() const
    {
        return _lateral->Direcao();
    }
private:
    TSuperficieCilindrica* _lateral = nullptr;
    TSuperficieCircular* _base = nullptr;
    TSuperficieCircular* _topo = nullptr;
};

// ------------------------------------------------------------------------------------------------

class TSuperficieConica : public IEntidade3D
{
public:
    TSuperficieConica(const TPonto3D& pCentroBase, double raioBase, double altura, const TVetor3D& direcao)
        : _c(pCentroBase), _r(raioBase), _h(altura), _d(direcao.Normalizado())
    {
        _v = _c + _d * _h;
        _cosTetaPow2 = (_h * _h) / (_h * _h + _r * _r);
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

    TMaterial Material(const TRaio3D&) const override
    {
        return _material;
    }
    void Material(const TMaterial& material)
    {
        _material = material;
    }

    TCoordenadasUV CoordenadasUV(const TPonto3D& p, const TRaio3D& raio) const override
    {
        return { 0.0, 0.0 };
    }

    // nos dois metodos abaixo, as notas de aula estavam bem diferentes
    // das fotos da lousa. inicialmente tentei implementar pelas fotos,
    // mas ocorreram diversos problemas, entao fui pelas notas, que
    // pareceu funcionar com alguns ajustes. depois avaliar o que errei
    // ao fazer a versao da lousa (p.s: aparentemente era bem mais
    // pesada, pelo uso ostensivo de matrizes)
    TVetor3D Normal(const TPonto3D& p, const TRaio3D& r) const override
    {
        const TVetor3D w = p - _v;

        return TVetor3D { _d * _d.Dot(w) - w * _cosTetaPow2 }.Normalizado();
    }
    std::vector<double> Intersecoes(const TRaio3D& raio) const override
    {
        std::vector<double> intersecoes;

        const TVetor3D& dr = raio.Direcao();
        const TVetor3D& dc = _d;

        const TVetor3D v = _v - raio.Origem();
        const double dn = dr.Dot(dc);
        const double dd = dr.Dot(dr);
        const double vn = v.Dot(dc);
        const double vd = v.Dot(dr);
        const double vv = v.Dot(v);

        const double a = dn * dn - dd * _cosTetaPow2;
        const double b = vd * _cosTetaPow2 - vn * dn;
        const double c = vn * vn - vv * _cosTetaPow2;

        // 1) ver aqui como escrever da forma tradicional, ex: 4ac,
        //    para poder passar para o solver generico depois
        // 2) talvez fosse bom avaliar "a" como fabs(a) < 1e-12
        const double delta = b * b - a * c;
        if (a == 0.0 || delta < 0.0)
        {
            return intersecoes;
        }

        // ver aqui como escrever da forma tradicional, ex: / 2.0a
        // para poder passar para o solver generico depois
        const double t1 = (-b - sqrt(delta)) / (a);
        const double t2 = (-b + sqrt(delta)) / (a);
        const TPonto3D p1 = raio.Ponto(t1);
        const TPonto3D p2 = raio.Ponto(t2);
        const double proj1 = dc.Dot(_v - p1);
        const double proj2 = dc.Dot(_v - p2);
        const bool t1Valido = FuncoesGerais::EstaEm(proj1, 0.0, _h);
        const bool t2Valido = FuncoesGerais::EstaEm(proj2, 0.0, _h);

        if (t1Valido)
        {
            intersecoes.push_back(t1);
        }
        if (t2Valido)
        {
            intersecoes.push_back(t2);
        }

        // nao entendi pq foi necessario ordenar, acho que ja devia vir
        // na ordem certa, mas pelo visto resolveu, entao devia estar
        // vindo um t maior no comeco, avaliar pq
        std::sort(intersecoes.begin(), intersecoes.end());

        return intersecoes;
    }

    void Transforma(const TMatriz<double>& matriz4x4) override
    {
        // ToDo
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
    double _cosTetaPow2 = 0.0;
};

// ------------------------------------------------------------------------------------------------

class TCone : public TEntidadeComposta
{
public:
    TCone(const TPonto3D& pCentroBase, double raioBase, double altura, const TVetor3D& direcao)
    {
        Insere(TSuperficieConica { pCentroBase, raioBase, altura, direcao });
        Insere(TSuperficieCircular { pCentroBase, direcao, raioBase });

        _lateral = static_cast<TSuperficieConica*>(_entidades[0].get());
        _base = static_cast<TSuperficieCircular*>(_entidades[1].get());
    }

    void Material(const TMaterial& material)
    {
        _lateral->Material(material);
        _base->Material(material);
    }

    const TPonto3D& CentroBase() const
    {
        return _lateral->CentroBase();
    }
    double RaioBase() const
    {
        return _lateral->RaioBase();
    }
    double Altura() const
    {
        return _lateral->Altura();
    }
    const TVetor3D& Direcao() const
    {
        return _lateral->Direcao();
    }
private:
    TSuperficieConica* _lateral = nullptr;
    TSuperficieCircular* _base = nullptr;
};

// ------------------------------------------------------------------------------------------------

class TLeitorMalha3D
{
private:
    struct TVertice
    {
        int v  = -1; // indice da posicao geometrica do vertice
        int vt = -1; // indice da coordenada da textura do vertice
        int vn = -1; // indice da normal do vertice
    };

    struct TFace
    {
        std::vector<TVertice> vertices;
        TMaterial* material = nullptr;
    };

public:
    TLeitorMalha3D(const std::string& obj) : _caminhoObj(obj), _isObj(obj) {}

    static void Carrega(const std::string& obj, TMalha3D& malha)
    {
        const bool bkpDeveRotular = malha.DeveRotularTriangulosInseridosAutomaticamente();
        malha.DeveRotularTriangulosInseridosAutomaticamente(true);

        TLeitorMalha3D leitorMalha(obj);
        leitorMalha.Parse();
        leitorMalha.Popula(malha);

        malha.DeveRotularTriangulosInseridosAutomaticamente(bkpDeveRotular);
    }

    bool Parse()
    {
        if (_parseado)
        {
            return _ok;
        }

        if (_isObj.is_open() && PreProcessaOBJ())
        {
            _isMtl.open(_caminhoMtl);
        }
        else
        {
            _ok = false;
        }

        if (_isMtl.is_open())
        {
            ProcessaMTL();
        }
        else
        {
            _ok = false;
        }

        if (_ok)
        {
            _isObj.seekg(0, std::ios_base::beg);
            ProcessaOBJ();
        }

        _parseado = true;
        
        return _ok;
    }

    bool Ok() const
    {
        return _ok;
    }

    void Popula(TMalha3D& malha) const
    {
        // Eh uma premissa de 99% dos exporters, incluindo obj, que as
        // faces declaradas sao poligonos planos convexos. Se nao for,
        // nem o proprio OpenGL aceita. Portanto, para qualquer face com N
        // vertices (N >= 3), podemos triangularizar da seguinte forma:
        // (v0, v1, v2)
        // (v0, v2, v3)
        // (v0, v3, v4)
        // ...
        // (v0, v(N-2), v(N-1))
        // , onde temos um total de N - 2 triangulos
        // 
        // Por exemplo, considerando N = 5, temos NT = 5 - 2 = 3, a saber
        // T1 = (v0, v1, v2)
        // T2 = (v0, v2, v3)
        // T3 = (v0, v3, v4)

        for (const TFace& face : _faces)
        {
            const std::vector<TVertice>& verticesFace = face.vertices;
            const auto nVertices = static_cast<int>(verticesFace.size());
            if (nVertices < 3)
            {
                continue;
            }

            const int nt = nVertices - 2;
            for (int i = 0; i < nt; i++)
            {
                const TVertice& v1 = verticesFace[0];
                const TVertice& v2 = verticesFace[1 + i];
                const TVertice& v3 = verticesFace[2 + i];

                const TPonto3D& p1 = _verticesGeometricos[v1.v];
                const TPonto3D& p2 = _verticesGeometricos[v2.v];
                const TPonto3D& p3 = _verticesGeometricos[v3.v];

                TSuperficieTriangular t { p1, p2, p3 };

                if (v1.vt != -1)
                {
                    const TPonto3D& uv1 = _coordsTexturaVertices[v1.vt];
                    t.CoordenadasUV1({ uv1.X(), uv1.Y() });
                }

                if (v2.vt != -1)
                {
                    const TPonto3D& uv2 = _coordsTexturaVertices[v2.vt];
                    t.CoordenadasUV2({ uv2.X(), uv2.Y() });
                }

                if (v3.vt != -1)
                {
                    const TPonto3D& uv3 = _coordsTexturaVertices[v3.vt];
                    t.CoordenadasUV3({ uv3.X(), uv3.Y() });
                }

                if (face.material != nullptr)
                {
                    t.Material(*face.material);
                }

                malha.Adiciona(t);
            }
        }
    }

private:
    bool PreProcessaOBJ()
    {
        bool achouMtlAssociado = false;

        std::string linha;
        while (std::getline(_isObj, linha))
        {
            NormalizaLinha(linha);

            if (IgnoraLinha(linha))
            {
                // NoOp
            }
            else if (IniciaComPalavraChave(linha, "mtllib"))
            {
                _caminhoMtl = CaminhoArquivoNoMesmoDiretorio(LeString(linha));
                achouMtlAssociado = true;
                break;
            }
        }

        return achouMtlAssociado;
    }

    void ProcessaMTL()
    {
        std::string linha;
        while (std::getline(_isMtl, linha))
        {
            NormalizaLinha(linha);
            ProcessaLinhaMTL(linha);
        }
    }

    void ProcessaOBJ()
    {
        std::string linha;
        while (std::getline(_isObj, linha))
        {
            NormalizaLinha(linha);
            ProcessaLinhaOBJ(linha);
        }
    }

    void NormalizaLinha(std::string& linha) const
    {
        size_t inicio = linha.find_first_not_of(" \t");
        if (inicio != std::string::npos)
        {
            linha = linha.substr(inicio);
        }
    }

    void ProcessaLinhaMTL(const std::string& linha)
    {
        if (IgnoraLinha(linha)) // Linha em branco ou comentario
        {
            // NoOp
        }
        else if (IniciaComPalavraChave(linha, "newmtl")) // Nome Material
        {
            _materialCorrente = LeString(linha);
        }
        else if (IniciaComPalavraChave(linha, "Ns")) // Brilho
        {
            TMaterial& materialCorrente = MaterialCorrente();

            const double m = LeEscalar(linha);
            materialCorrente.M(m);
        }
        else if (IniciaComPalavraChave(linha, "Ka")) // Fator Luz Ambiente
        {
            TMaterial& materialCorrente = MaterialCorrente();

            const TVetor3D ka = LeVetor3D(linha);
            materialCorrente.KaR(ka.X());
            materialCorrente.KaG(ka.Y());
            materialCorrente.KaB(ka.Z());
        }
        else if (IniciaComPalavraChave(linha, "Kd")) // Fator Luz Difusa
        {
            TMaterial& materialCorrente = MaterialCorrente();

            const TVetor3D kd = LeVetor3D(linha);
            materialCorrente.KdR(kd.X());
            materialCorrente.KdG(kd.Y());
            materialCorrente.KdB(kd.Z());
        }
        else if (IniciaComPalavraChave(linha, "Ks")) // Fator Luz Especular
        {
            TMaterial& materialCorrente = MaterialCorrente();

            const TVetor3D ke = LeVetor3D(linha);
            materialCorrente.KeR(ke.X());
            materialCorrente.KeG(ke.Y());
            materialCorrente.KeB(ke.Z());
        }
        else if (IniciaComPalavraChave(linha, "Ke")) // Fator Luz Emissao
        {
            // const TVetor3D k = LeVetor3D(linha);
        }
        else if (IniciaComPalavraChave(linha, "Ni")) // Densidade Optica
        {
            // const double ni = LeEscalar(linha);
        }
        else if (IniciaComPalavraChave(linha, "d")) // Transparencia
        {
            // const double d = LeEscalar(linha);
        }
        else if (IniciaComPalavraChave(linha, "illum")) // Modelo Iluminacao
        {
            // const int illum = LeInteiro(linha);
        }
        else if (IniciaComPalavraChave(linha, "map_Kd")) // Mapa Textura
        {
            TMaterial& materialCorrente = MaterialCorrente();

            const std::string mapKd = LeString(linha);
            materialCorrente.CarregaTextura(CaminhoArquivoNoMesmoDiretorio(mapKd));
        }
    }

    void ProcessaLinhaOBJ(const std::string& linha)
    {
        if (IgnoraLinha(linha)) // Linha em branco ou comentario
        {
            // NoOp
        }
        else if (IniciaComPalavraChave(linha, "mtllib")) // Nome Arquivo Materiais
        {
            // NoOp, ja foi tratado antes para carregar o MTL
        }
        else if (IniciaComPalavraChave(linha, "usemtl")) // Nome Material
        {
            _materialCorrente = LeString(linha);
        }
        else if (IniciaComPalavraChave(linha, "v")) // Vertice Geometrico
        {
            _verticesGeometricos.push_back(LeVetor3D(linha));
        }
        else if (IniciaComPalavraChave(linha, "vt")) // Coordenadas Textura Vertice
        {
            _coordsTexturaVertices.push_back(LeVetor2D(linha));
        }
        else if (IniciaComPalavraChave(linha, "vn")) // Normal Vertice
        {
            _normaisVertices.push_back(LeVetor3D(linha));
        }
        else if (IniciaComPalavraChave(linha, "f")) // Face
        {
            ProcessaFace(linha);
        }
    }

    void ProcessaFace(const std::string& linha)
    {
        auto LeVertice = [](const std::string& vertice)
        {
            std::istringstream iss(vertice);
            std::string parte;

            TVertice v;

            if (std::getline(iss, parte, '/'))
            {
                v.v = std::stoi(parte) - 1;
            }

            if (std::getline(iss, parte, '/') && !parte.empty())
            {
                v.vt = std::stoi(parte) - 1;
            }

            if (std::getline(iss, parte, '/') && !parte.empty())
            {
                v.vn = std::stoi(parte) - 1;
            }

            return v;
        };

        std::istringstream iss(linha);

        std::string token;
        iss >> token;

        std::vector<TVertice> verticesFace;

        std::string strVertice;
        while (iss >> strVertice)
        {
            const TVertice vertice = LeVertice(strVertice);
            verticesFace.push_back(vertice);
        }
        
        TFace face;
        face.vertices = verticesFace;
        face.material = &MaterialCorrente();

        _faces.push_back(face);
    }

    std::string LeString(const std::string& linha) const
    {
        std::string token;
        std::string str;

        std::istringstream iss(linha);
        iss >> token >> str;

        return str;
    }

    int LeInteiro(const std::string& linha) const
    {
        std::string token;
        int n;

        std::istringstream iss(linha);
        iss >> token >> n;

        return n;
    }

    double LeEscalar(const std::string& linha) const
    {
        std::string token;
        double n;

        std::istringstream iss(linha);
        iss >> token >> n;

        return n;
    }

    TVetor3D LeVetor2D(const std::string& linha) const
    {
        std::string token;
        double x, y;

        std::istringstream iss(linha);
        iss >> token >> x >> y;

        return { x, y, 0.0 };
    }

    TVetor3D LeVetor3D(const std::string& linha) const
    {
        std::string token;
        double x, y, z;

        std::istringstream iss(linha);
        iss >> token >> x >> y >> z;

        return { x, y, z };
    }

    bool IgnoraLinha(const std::string& linha)
    {
        return linha.empty() || linha[0] == '#';
    }

    bool IniciaComPalavraChave(const std::string& linha, const std::string& token)
    {
        return linha.rfind(token + " ", 0) == 0;
    }

    TMaterial& MaterialCorrente()
    {
        return _materiais[_materialCorrente];
    }

    std::string CaminhoArquivoNoMesmoDiretorio(const std::string& arq) const
    {
        std::string caminhoArq;

        const size_t posBarra = _caminhoObj.find_last_of("/\\");
        if (posBarra == std::string::npos)
        {
            caminhoArq = arq;
        }
        else
        {
            const std::string diretorio = _caminhoObj.substr(0, posBarra + 1);
            caminhoArq = diretorio + arq;
        }

        return caminhoArq;
    }

    std::string _caminhoObj;
    std::string _caminhoMtl;
    std::ifstream _isObj;
    std::ifstream _isMtl;
    bool _parseado = false;
    bool _ok = true;

    std::unordered_map<std::string, TMaterial> _materiais;
    std::string _materialCorrente;

    std::vector<TPonto3D> _verticesGeometricos;
    std::vector<TPonto3D> _coordsTexturaVertices;
    std::vector<TVetor3D> _normaisVertices;
    std::vector<TFace> _faces;
};

// ------------------------------------------------------------------------------------------------

class TLeitorCena3D
{
public:
    TLeitorCena3D(const std::string& caminho) : _is(caminho) {}

private:
    std::ifstream _is; // .scn
};

// ------------------------------------------------------------------------------------------------

class TCamera
{
public:
    TCamera() = default;

    TCamera& Janela(double w, double h)
    {
        _wJanela = w;
        _hJanela = h;

        return *this;
    }
    TCamera& Viewport(uint16_t w, uint16_t h)
    {
        _wCanvas = w;
        _hCanvas = h;

        return *this;
    }
    TCamera& DistanciaFocal(double d)
    {
        _d = d;

        return *this;
    }
    TCamera& OlhoObservador(const TPonto3D& pEye)
    {
        _pEye = pEye;

        return *this;
    }
    TCamera& Visada(const TPonto3D& at)
    {
        _atPoint = at;

        return *this;
    }
    TCamera& Cima(const TVetor3D& up)
    {
        _lookUp = up;

        return *this;
    }

    TCamera& Init()
    {
        _w = TVetor3D(_pEye - _atPoint).Normalizado();
        _u = _lookUp.Vetorial(_w).Normalizado();
        _v = _w.Vetorial(_u);

        _dx = _wJanela / _wCanvas;
        _dy = _hJanela / _hCanvas;

        return *this;
    }

    double Largura() const { return _wJanela; }
    double Altura() const { return _hJanela; }
    uint16_t LarguraCanvas() const { return _wCanvas; }
    uint16_t AlturaCanvas() const { return _hCanvas; }

    TRaio3D RaioMundo(uint16_t x, uint16_t y) const
    {
        const double px = -0.5 * _wJanela + 0.5 * _dx + x * _dx;
        const double py =  0.5 * _hJanela - 0.5 * _dy - y * _dy;
        const double pz = -_d;

        const TPonto3D pCamera { px, py, pz };
        const TPonto3D pOrigemCamera { 0.0, 0.0, 0.0 }; // p0
        const TVetor3D direcaoRaioCamera = FuncoesGeometricas::Versor(pOrigemCamera, pCamera);

        const TVetor3D u_ = _u * direcaoRaioCamera.X();
        const TVetor3D v_ = _v * direcaoRaioCamera.Y();
        const TVetor3D w_ = _w * direcaoRaioCamera.Z();

        const TPonto3D& pOrigemMundo = _pEye;
        const TVetor3D direcaoRaioMundo = u_ + v_ + w_;

        return TRaio3D { pOrigemMundo, direcaoRaioMundo };
    }

private:
    double _wJanela = 0.0; // define o campo de visao horizontal (fovX)
    double _hJanela = 0.0; // define o campo de visao vertical   (fovY)
    uint16_t _wCanvas = 0; // largura da viewport
    uint16_t _hCanvas = 0; // altura  da viewport

    TPonto3D _pEye;    // olho do observador / posicao da camera
    TPonto3D _atPoint; // direcionamento de visada
    TVetor3D _lookUp;  // orientacao da camera em torno do eixo de visada
    double _d = 0.0;   // distancia focal

    TVetor3D _w; // z da camera, sentido positivo para tras dela
    TVetor3D _u; // x da camera, calculado pelo z dela e o upVector fornecido
    TVetor3D _v; // y da camera, pela regra da mao direita

    double _dx = 0.0;
    double _dy = 0.0;
};

// ------------------------------------------------------------------------------------------------

struct THit
{
    IEntidade3D* entidade;
    double interseccao;
};

// ------------------------------------------------------------------------------------------------

class TCena3D
{
public:
    TCena3D() = default;
    TCena3D(const TPonto3D& origem, const TCamera& camera) : _camera(camera), _p0(origem) {}

    const TCamera& Camera() const { return _camera; }
    void Camera(const TCamera& camera) { _camera = camera; }

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

    void Renderiza(IDispositivoSaida& arq)
    {
        if (_podeRenderizarMT)
        {
            RenderizaMultiThread(arq);
        }
        else
        {
            RenderizaSingleThread(arq);
        }
    }

    std::string Pick(uint16_t x, uint16_t y)
    {
        const TRaio3D raio = _camera.RaioMundo(x, y);
        const THit hit = ColisaoMaisProxima(raio);

        if (hit.entidade == nullptr)
        {
            return "NENHUM";
        }

        if (auto entidadeComposta = dynamic_cast<const TEntidadeComposta*>(hit.entidade))
        {
            return entidadeComposta->Rotulo(raio);
        }

        return hit.entidade->Rotulo();
    }

    // Nao sei se "colisao" eh uma boa traducao/adaptacao para hit, mas enfim
    // pensei em deixar em ingles, NearestHit ou ClosestHit, ja misturei mesmo
    THit ColisaoMaisProxima(const TRaio3D& raio) const
    {
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

        return { entidadeMaisProximaObservador, intersecaoMaisProximaObservador };
    }

    void Log(const std::string& msg) const
    {
        if (_arqLog != nullptr)
        {
            _arqLog->Anexa(msg);
        }
    }

    void PodeRenderizarMultiThread(bool mt)
    {
        _podeRenderizarMT = mt;
    }

private:
    void RenderizaSingleThread(IDispositivoSaida& arq)
    {
        if (auto arqLog = dynamic_cast<TArquivoLOG*>(&arq))
        {
            _arqLog = arqLog;
        }

        const uint16_t nLinhas = _camera.AlturaCanvas();
        const uint16_t nColunas = _camera.LarguraCanvas();

        for (int l = 0; l < nLinhas; l++)
        {
            for (int c = 0; c < nColunas; c++)
            {
                arq.Anexa(Cor(Camera().RaioMundo(c, l)));
            }
        }

        _arqLog = nullptr;
        arq.Flush();
    }

    void RenderizaMultiThread(IDispositivoSaida& arq)
    {
        struct RenderJob
        {
            RenderJob(IDispositivoSaida& arq, TCena3D& c, int y0, int y1)
                : arq(arq), c(c), y0(y0), y1(y1) {};

            void operator()()
            {
                for (int y = y0; y < y1; y++)
                    for (int x = 0; x < W(); x++)
                        arq.Anexa(c.Cor(c.Camera().RaioMundo(x, y)));
            }

            uint16_t W() const
            {
                return c.Camera().LarguraCanvas();
            }

            IDispositivoSaida& arq;
            TCena3D& c;
            int y0, y1;
        };

        const uint16_t h = _camera.AlturaCanvas();
        const int n = std::thread::hardware_concurrency();
        const int linhas = h / n;

        std::vector<std::thread> threads;
        for (int i = 0; i < n; i++)
        {
            int y0 = i * linhas;
            int y1 = (i == n - 1) ? h : y0 + linhas;
            threads.emplace_back(RenderJob { arq, *this, y0, y1 });
        }

        for (auto& t : threads) t.join();

        // const uint16_t nLinhas = _camera.AlturaCanvas();
        // const uint16_t nColunas = _camera.LarguraCanvas();

        // const double z = _camera.Centro().Z();
        // for (int l = 0; l < nLinhas; l++)
        // {
        //     const double y = _camera.Y(l);

        //     for (int c = 0; c < nColunas; c++)
        //     {
        //         const double x = _camera.X(c);

        //         arq.Anexa({ 255u, 0u, 0u });
        //     }
        // }

        // arq.Flush();
    }

    TCor Cor(const TRaio3D& raio) const
    {
        const THit hit = ColisaoMaisProxima(raio);
        const IEntidade3D* entidadeMaisProximaObservador = hit.entidade;
        const double intersecaoMaisProximaObservador = hit.interseccao;

        TCor pixel = _bgColor;

        if (entidadeMaisProximaObservador != nullptr)
        {
            pixel = Cor(*entidadeMaisProximaObservador, raio, intersecaoMaisProximaObservador);
        }

        return pixel;
    }

    // Segundo o ChatGPT, essa minha funcao aqui seria meu "shader", escrito "na mao"
    // Ouco muito falar em shader, mas nao entendo bem o que eh. Grosso modo, ele
    // disse que eh a funcao que calcula a cor final de um pixel (renderiza um pixel?),
    // entao no meu caso seria essa junto com a de cima (mas principalmente essa)
    TCor Cor(const IEntidade3D& entidade, const TRaio3D& raio, double ti) const
    {
        const TPonto3D pi = raio.Ponto(ti);
        const TCoordenadasUV uv = entidade.CoordenadasUV(pi, raio);

        const TMaterial& material = entidade.Material(raio);
        const double kdR = material.KdR(uv.u, uv.v);
        const double kdG = material.KdG(uv.u, uv.v);
        const double kdB = material.KdB(uv.u, uv.v);
        const TPonto3D kd { kdR, kdG, kdB };
        const TPonto3D ke { material.KeR(), material.KeG(), material.KeB() };
        const TPonto3D ka { material.KaR(), material.KaG(), material.KaB() };

        const TVetor3D iAmb { IambR(), IambG(), IambB() };
        const auto ia = iAmb.Arroba(ka);

        const TVetor3D n = entidade.Normal(pi, raio);
        const TVetor3D v = raio.Direcao() * -1.0;

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
    bool _podeRenderizarMT = false;

    TCamera _camera;
    TPonto3D _p0; // olho do pintor (origem)

    TCor _bgColor;
    double _iAmbR = 0.0;
    double _iAmbG = 0.0;
    double _iAmbB = 0.0;

    std::vector<std::unique_ptr<IEntidade3D>> _entidades;
    std::vector<std::unique_ptr<IFonteLuminosa>> _fontes;
};

// ------------------------------------------------------------------------------------------------

TMaterial FabricaMaterialHomogeneo(const TVetor3D& k, double m)
{
    TMaterial material;

    material.KdR(k.X());
    material.KdG(k.Y());
    material.KdB(k.Z());
    material.KeR(k.X());
    material.KeG(k.Y());
    material.KeB(k.Z());
    material.KaR(k.X());
    material.KaG(k.Y());
    material.KaB(k.Z());
    material.M(m);

    return material;
}

// ------------------------------------------------------------------------------------------------

TPlano FabricaChao()
{
    TMaterial material;
    material.KdR(0.0);
    material.KdG(0.0);
    material.KdB(0.0);
    material.KeR(0.0);
    material.KeG(0.0);
    material.KeB(0.0);
    material.KaR(0.0);
    material.KaG(0.0);
    material.KaB(0.0);
    material.M(1.0);
    material.CarregaTextura(ResTbl::TEX_MADEIRA);
    material.Textura()->K(30.0);

    TPlano planoChao { { 0.0, -150.0, 0.0 }, { 0.0, 1.0, 0.0 } };
    planoChao.Rotulo("PLANO_CHAO");
    planoChao.Material(material);

    return planoChao;
}

// ------------------------------------------------------------------------------------------------

TPlano FabricaParedeLateralDireita()
{
    TPlano planoParedeDireita { { 200.0, -150.0, 0.0 }, { -1.0, 0.0, 0.0 } };
    planoParedeDireita.Rotulo("PLANO_PAREDE_LATERAL_DIR");
    planoParedeDireita.Material(FabricaMaterialHomogeneo({ 0.686, 0.933, 0.933 }, 1.0));

    return planoParedeDireita;
}

// ------------------------------------------------------------------------------------------------

TPlano FabricaParedeFrontal()
{
    TPlano planoParedeFrontal { { 200.0, -150.0, -400.0 }, { 0.0, 0.0, 1.0 } };
    planoParedeFrontal.Rotulo("PLANO_PAREDE_FRONTAL");
    planoParedeFrontal.Material(FabricaMaterialHomogeneo({ 0.686, 0.933, 0.933 }, 1.0));

    return planoParedeFrontal;
}

// ------------------------------------------------------------------------------------------------

TPlano FabricaParedeLateralEsquerda()
{
    TPlano planoParedeEsquerda { { -200.0, -150.0, 0.0 }, { 1.0, 0.0, 0.0 } };
    planoParedeEsquerda.Rotulo("PLANO_PAREDE_LATERAL_ESQ");
    planoParedeEsquerda.Material(FabricaMaterialHomogeneo({ 0.686, 0.933, 0.933 }, 1.0));

    return planoParedeEsquerda;
}

// ------------------------------------------------------------------------------------------------

TPlano FabricaTeto()
{
    TPlano planoTeto { { 0.0, 150.0, 0.0 }, { 0.0, -1.0, 0.0 } };
    planoTeto.Rotulo("PLANO_TETO");
    planoTeto.Material(FabricaMaterialHomogeneo({ 0.933, 0.933, 0.933 }, 1.0));

    return planoTeto;
}

// ------------------------------------------------------------------------------------------------

TCilindro FabricaCilindro()
{
    const TPonto3D cBaseCilindro { 0.0, -150.0, -200.0 };
    const double rBaseCilindro = 5.0;
    const double hCilindro = 90.0;
    const TVetor3D dCilindro { 0.0, 1.0, 0.0 };

    TCilindro cilindro { cBaseCilindro, rBaseCilindro, hCilindro, dCilindro };
    cilindro.Rotulo("CILINDRO_1");
    cilindro.Material(FabricaMaterialHomogeneo({ 0.824, 0.706, 0.549 }, 10.0));

    return cilindro;
}

// ------------------------------------------------------------------------------------------------

TCone FabricaCone()
{
    const TPonto3D cBaseCone { 0.0, -60.0, -200.0 };
    const double rBaseCone = 90.0;
    const double hCone = 150.0;
    const TVetor3D& dCone = { 0.0, 1.0, 0.0 };

    TCone cone { cBaseCone, rBaseCone, hCone, dCone };
    cone.Rotulo("CONE_1");
    cone.Material(FabricaMaterialHomogeneo({ 0.0, 1.0, 0.498 }, 10.0));

    return cone;
}

// ------------------------------------------------------------------------------------------------

TCubo FabricaCubo()
{
    TCubo cubo { { 0.0, -150.0, -165.0 }, 40.0 };
    cubo.Rotulo("CUBO_1");
    cubo.Material(FabricaMaterialHomogeneo({ 1.0, 0.078, 0.576 }, 10.0));

    return cubo;
}

// ------------------------------------------------------------------------------------------------

TMalha3D FabricaMalha()
{
    TMalha3D malha;
    malha.Rotulo("MALHA_1");
    TLeitorMalha3D::Carrega(ResTbl::OBJ_SPYRO, malha);

    return malha;
}

// ------------------------------------------------------------------------------------------------

TEsfera FabricaEsfera()
{
    TEsfera esfera { { 0.0, 95.0, -200.0 }, 5.0 };
    esfera.Rotulo("ESFERA_1");
    esfera.Material(FabricaMaterialHomogeneo({ 0.854, 0.647, 0.125 }, 10.0));

    return esfera;
}

// ------------------------------------------------------------------------------------------------

TFontePontual FabricaFontePontual()
{
    return { { -100.0, 140.0, -20.0 }, { 0.7, 0.7, 0.7 } };
}

// ------------------------------------------------------------------------------------------------

TCamera FabricaCamera()
{
    // Constroi uma camera identidade (espaco de camera = espaco de mundo),
    // equivalente a que viemos usando desde o trabalho 1

    return TCamera()
            .Janela(60.0, 60.0)
            .DistanciaFocal(30.0)
            .Viewport(500u, 500u)
            .OlhoObservador({ 0.0, 0.0, 0.0 })
            .Visada({ 0.0, 0.0, -1.0 })
            .Cima({ 0.0, 1.0, 0.0 })
            .Init();
}

// ------------------------------------------------------------------------------------------------

TCena3D FabricaCena1()
{
    const TPonto3D p0 { 0.0, 0.0, 0.0 };

    TCena3D cena { p0, FabricaCamera() };
    cena.BgColor({ 100u, 100u, 100u });
    cena.IambR(0.3);
    cena.IambG(0.3);
    cena.IambB(0.3);

    cena.Insere(FabricaChao());
    cena.Insere(FabricaParedeLateralDireita());
    cena.Insere(FabricaParedeFrontal());
    cena.Insere(FabricaParedeLateralEsquerda());
    cena.Insere(FabricaTeto());
    cena.Insere(FabricaCilindro());
    cena.Insere(FabricaCone());
    cena.Insere(FabricaCubo());
    cena.Insere(FabricaEsfera());
    cena.Insere(FabricaFontePontual());

    cena.PodeRenderizarMultiThread(false);

    return cena;
}

// ------------------------------------------------------------------------------------------------

TCena3D FabricaCena2()
{
    const TPonto3D p0 { 0.0, 0.0, 0.0 };

    TCena3D cena { p0, FabricaCamera() };
    cena.BgColor({ 100u, 100u, 100u });
    cena.IambR(0.3);
    cena.IambG(0.3);
    cena.IambB(0.3);

    cena.Insere(FabricaChao());
    cena.Insere(FabricaParedeLateralDireita());
    cena.Insere(FabricaParedeFrontal());
    cena.Insere(FabricaParedeLateralEsquerda());
    cena.Insere(FabricaTeto());
    cena.Insere(FabricaMalha());
    cena.Insere(FabricaFontePontual());

    cena.PodeRenderizarMultiThread(false);

    return cena;
}

// ------------------------------------------------------------------------------------------------

TCena3D FabricaCena3()
{
    const TPonto3D p0 { 0.0, 0.0, 0.0 };

    TCena3D cena { p0, FabricaCamera().Janela(30.0, 30.0).Viewport(250u, 250u).Init() };
    cena.BgColor({ 100u, 100u, 100u });
    cena.IambR(0.3);
    cena.IambG(0.3);
    cena.IambB(0.3);

    cena.Insere(FabricaChao());
    cena.Insere(FabricaParedeLateralDireita());
    cena.Insere(FabricaParedeFrontal());
    cena.Insere(FabricaParedeLateralEsquerda());
    cena.Insere(FabricaTeto());
    cena.Insere(FabricaMalha());
    cena.Insere(TFontePontual { { -4.0, 4.0, 0.0 }, { 0.7, 0.7, 0.7 } });
    // cena.Insere(FabricaFontePontual());

    cena.PodeRenderizarMultiThread(false);

    return cena;
}

// ------------------------------------------------------------------------------------------------

TCena3D FabricaCena4()
{
    const TPonto3D p0 { 0.0, 0.0, 0.0 };

    TCena3D cena { p0, FabricaCamera().Janela(15.0, 15.0).Viewport(125u, 125u).Init() };
    cena.BgColor({ 100u, 100u, 100u });
    cena.IambR(0.3);
    cena.IambG(0.3);
    cena.IambB(0.3);

    // cena.Insere(FabricaChao());
    // cena.Insere(FabricaParedeLateralDireita());
    // cena.Insere(FabricaParedeFrontal());
    // cena.Insere(FabricaParedeLateralEsquerda());
    // cena.Insere(FabricaTeto());
    cena.Insere(FabricaMalha());
    cena.Insere(FabricaFontePontual());

    cena.PodeRenderizarMultiThread(false);

    return cena;
}

// ------------------------------------------------------------------------------------------------

typedef TCena3D(*PtrFabricaCena)();

std::unordered_map<std::string, PtrFabricaCena> FabricasCenasPreDefinidas = {
    { "NATAL", FabricaCena1 },
    { "SPYRO", FabricaCena2 },
    { "TESTE", FabricaCena3 },
    { "MENOR", FabricaCena4 }
};

// ------------------------------------------------------------------------------------------------

TCena3D FabricaCena()
{
    return FabricasCenasPreDefinidas["TESTE"]();
}

// ------------------------------------------------------------------------------------------------

std::unique_ptr<IArquivoSaida> FabricaArquivo(
    const std::string& nome,
    EFormatoImagem formato,
    const TCena3D& cena
)
{
    const TCamera camera = cena.Camera();

    return FuncoesGerais::FabricaArquivo(
        formato, nome, camera.LarguraCanvas(), camera.AlturaCanvas()
    );
}

// ------------------------------------------------------------------------------------------------

bool RenderizaImagem()
{
    TCena3D cena = FabricaCena();
    const std::unique_ptr<IArquivoSaida> arq = FabricaArquivo("a", EFormatoImagem::BMP, cena);

    const bool erro = !arq->Aberto();
    if (!erro)
    {
        cena.Renderiza(*arq);
    }

    return erro;
}

// ------------------------------------------------------------------------------------------------

class Globals
{
public:
    Globals(const Globals&) = delete;
    Globals(Globals&&) = delete;
    Globals& operator=(const Globals&) = delete;

    ~Globals()
    {
        delete cena;
    }

    static Globals& Instancia()
    {
        if (_instancia == nullptr)
        {
            _instancia = new Globals;
        }

        return *_instancia;
    }

    TCena3D& Cena()
    {
        return *cena;
    }

    std::wstring PickObject(int x, int y)
    {
        const auto _x = static_cast<std::uint16_t>(x);
        const auto _y = static_cast<std::uint16_t>(y);
        const std::string rotulo = cena->Pick(_x, _y);

        std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>> converter;
        const std::wstring pickedObj = converter.from_bytes(rotulo);

        return pickedObj;
    }

private:
    Globals()
    {
        cena = new TCena3D(FabricaCena());
    }

    TCena3D* cena = nullptr;

    static Globals* _instancia;
};

Globals* Globals::_instancia = nullptr;

// ------------------------------------------------------------------------------------------------

class TMainWindow
{
public:
    TMainWindow(
        HINSTANCE hInstance,
        PWSTR pCmdLine,
        int nCmdShow
    ) :
        _hInstance(hInstance),
        _pCmdLine(pCmdLine),
        _nCmdShow(nCmdShow),
        _clsName(TituloJanela()),
        _wndTitle(TituloJanela())
    {
    };

    static constexpr const wchar_t* TituloJanela()
    {
        return L"CG1";
    }
    static constexpr std::uint16_t XAncoraViewport()
    {
        return 160;
    }
    static constexpr std::uint16_t YAncoraViewport()
    {
        return 16;
    }

    void Executa()
    {
        if (MontaJanelaPrincipal())
        {
            ExibeJanelaPrincipal();
            ProcessaMensagens();
        }
    }

private:
    bool MontaJanelaPrincipal()
    {
        _wndCls = CriaClasse();
        RegisterClass(&_wndCls);

        _hWnd = CriaJanelaPrincipal();
        _hWndBtnRenderScene = CriaBotao(_hWnd, _IDC_BTN_RENDER_SCENE, L"Renderiza Cena", 16, 16, 128, 32);

        return _hWnd != nullptr && _hWndBtnRenderScene != nullptr;
    }

    WNDCLASS CriaClasse() const
    {
        WNDCLASS wc = {};
        wc.lpfnWndProc   = TMainWindow::WindowProc;
        wc.hInstance     = _hInstance;
        wc.lpszClassName = _clsName.c_str();

        return wc;
    }

    HWND CriaJanelaPrincipal() const
    {
        HWND hWnd = CreateWindowEx(
            0,                              // Optional window styles
            _clsName.c_str(),               // Window class
            _wndTitle.c_str(),              // Window text
            WS_OVERLAPPEDWINDOW,            // Window style

            // Size and position
            CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT, CW_USEDEFAULT,

            NULL,       // Parent window    
            NULL,       // Menu
            _hInstance, // Instance handle
            NULL        // Additional application data
        );

        return hWnd;
    }

    HWND CriaBotao(HWND parent, int id, const std::wstring& txt, int x, int y, int w, int h) const
    {
        HWND hWnd = CreateWindowEx( 
            0,                   // Optional window styles
            L"BUTTON",           // Predefined class; Unicode assumed
            txt.c_str(),         // Button text
            WS_TABSTOP | WS_VISIBLE | WS_CHILD | BS_DEFPUSHBUTTON, // Styles
            x,                   // x position
            y,                   // y position
            w,                   // Button width
            h,                   // Button height
            parent,              // Parent window
            (HMENU)(intptr_t)id, // Control ID for child windows
            (HINSTANCE)GetWindowLongPtr(parent, GWLP_HINSTANCE), 
            NULL);               // Pointer not needed.

        return hWnd;
    }

    void ExibeJanelaPrincipal()
    {
        ShowWindow(_hWnd, _nCmdShow);
    }

    static void EvBtnRenderScene(HWND hWndJanelaPrincipal)
    {
        rendering = true;

        InvalidateRect(hWndJanelaPrincipal, NULL, true);
        UpdateWindow(hWndJanelaPrincipal);
    }

    static void RenderScene(HDC hdc)
    {
        TCena3D& cena = Globals::Instancia().Cena();
        const std::uint16_t x = XAncoraViewport();
        const std::uint16_t y = YAncoraViewport();
        const std::uint16_t w = cena.Camera().LarguraCanvas();
        const std::uint16_t h = cena.Camera().AlturaCanvas();

        RECT r { x, y, x + w, y + h };
        FillRect(hdc, &r, (HBRUSH)GetStockObject(BLACK_BRUSH));

        if (rendering)
        {
            HCURSOR cursorAntigo = GetCursor();
            SetCursor(LoadCursor(NULL, IDC_WAIT));

            // Um teste interessante eh comparar o tempo de desenho propriamente dito
            // entre renderizar direto na viewport, ou no framebuffer como device
            // intermediario fazendo flush na viewport. O tempo total eh o mesmo, ja
            // que nao tem como fugir dos calculos pesados para cada pixel, mas eh
            // interessante comparar a velocidade de escrita dos pixels on demand
            // versus previamente armazenados no framebuffer
            // TFrameBuffer fb(w, h);
            TWin32Viewport vp(hdc, x, y, w, h);
            // fb.OutDevice(vp);
            // cena.Renderiza(fb);
            cena.Renderiza(vp);

            SetCursor(cursorAntigo);
        }

        rendering = false;
    }

    void ProcessaMensagens()
    {
        MSG msg = {};
        while (GetMessage(&msg, NULL, 0, 0))
        {
            TranslateMessage(&msg);
            DispatchMessage(&msg);
        }
    }

    static LRESULT CALLBACK WindowProc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
    {
        switch (uMsg)
        {
            case WM_DESTROY   : return EvQuit();
            case WM_PAINT     : return EvPaint(hWnd);
            case WM_SETCURSOR : if (EvSetCursor(lParam)) return true; break;
            case WM_MOUSEMOVE : EvMouseMove(hWnd, lParam); break;
            case WM_LBUTTONUP : EvMouseLeftButtonUp(hWnd, lParam); break;
            case WM_COMMAND   : return EvCommand(hWnd, wParam);
        }

        return DefWindowProc(hWnd, uMsg, wParam, lParam);
    }

    static bool EvQuit()
    {
        PostQuitMessage(0);

        return 0;
    }

    static bool EvPaint(HWND hWnd)
    {
        PAINTSTRUCT ps;
        HDC hdc = BeginPaint(hWnd, &ps);
        FillRect(hdc, &ps.rcPaint, (HBRUSH) (COLOR_WINDOW+1));
        RenderScene(hdc);
        EndPaint(hWnd, &ps);

        return 0;
    }

    static bool EvSetCursor(LPARAM lParam)
    {
        bool ok = false;
        if (LOWORD(lParam) == HTCLIENT)
        {
            SetCursor(LoadCursor(NULL, IDC_ARROW));
            ok = true;
        }

        return ok;
    }

    static bool EvMouseMove(HWND hWnd, LPARAM lParam)
    {
        const int xMousePos = GET_X_LPARAM(lParam);
        const int yMousePos = GET_Y_LPARAM(lParam);

        wchar_t buf[BUFLEN];
        GetClassName(hWnd, buf, BUFLEN);

        std::wstringstream ss;
        ss << buf << " ( X: " << xMousePos << ", Y: " << yMousePos << " ) - OBJ: " << pickedObject;
        auto str = ss.str();
        
        SetWindowText(hWnd, str.c_str());

        return true;
    }

    static bool EvMouseLeftButtonUp(HWND hWnd, LPARAM lParam)
    {
        int xMousePos = GET_X_LPARAM(lParam);
        int yMousePos = GET_Y_LPARAM(lParam);
        if (WithinViewportClientArea(xMousePos, yMousePos))
        {
            ConvertToViewportCoords(xMousePos, yMousePos);
            pickedObject = Globals::Instancia().PickObject(xMousePos, yMousePos);
        }
        else
        {
            pickedObject = L"NENHUM";
        }

        EvMouseMove(hWnd, lParam);

        return true;
    }

    static bool EvCommand(HWND hWnd, WPARAM wParam)
    {
        if (HIWORD(wParam) == BN_CLICKED)
        {
            if (LOWORD(wParam) == _IDC_BTN_RENDER_SCENE)
            {
                EvBtnRenderScene(hWnd);
            }
        }

        return 0;
    }

    static bool WithinViewportClientArea(std::uint16_t x, std::uint16_t y)
    {
        TCena3D& cena = Globals::Instancia().Cena();
        const auto xAncVp = XAncoraViewport();
        const auto yAncVp = YAncoraViewport();

        return x >= xAncVp && x < xAncVp + cena.Camera().LarguraCanvas() &&
               y >= yAncVp && y < yAncVp + cena.Camera().AlturaCanvas();
    }

    static void ConvertToViewportCoords(int& x, int& y)
    {
        if (WithinViewportClientArea(x, y))
        {
            x -= XAncoraViewport();
            y -= YAncoraViewport();
        }
    }

    std::wstring _clsName;
    std::wstring _wndTitle;
    HINSTANCE _hInstance;
    PWSTR _pCmdLine;
    int _nCmdShow;

    HWND _hWnd;
    WNDCLASS _wndCls;

    HWND _hWndBtnRenderScene; static constexpr int _IDC_BTN_RENDER_SCENE = 1;

    static bool rendering;
    static std::wstring pickedObject;
};

bool TMainWindow::rendering = false;
std::wstring TMainWindow::pickedObject = L"NENHUM";

// ------------------------------------------------------------------------------------------------

int WINAPI wWinMain(HINSTANCE hInstance, HINSTANCE /*hPrevInstance*/, PWSTR pCmdLine, int nCmdShow)
{
    Tracer::LimpaTrace();
    TMainWindow(hInstance, pCmdLine, nCmdShow).Executa();

    return 0;
}

// ------------------------------------------------------------------------------------------------
