#pragma once

#define INFO(m) \
    {\
    std::cerr << "\033[44m\033[1m[INFO]\033[49m\033[0m "; \
    std::cerr << m << std::endl; \
    }

#define WARNING(m) \
    {\
    std::cerr << std::endl; \
    std::cerr << "\033[43m\033[1m[WARNING]\033[49m\033[0m " << std::endl; \
    std::cerr << "\033[1mFILE\033[0m    : " << __FILE__ << std::endl; \
    std::cerr << "\033[1mFUNCTION\033[0m: " << __func__ << std::endl; \
    std::cerr << "\033[1mLINE\033[0m    : " << __LINE__ << std::endl; \
    std::cerr << "\033[1mMESSAGE\033[0m : " << m << std::endl; \
    std::cerr << std::endl; \
    }

#define ERROR(m) \
    {\
    std::cerr << std::endl; \
    std::cerr << "\033[41m\033[1m[ERROR]\033[49m\033[0m " << std::endl; \
    std::cerr << "\033[1mFILE\033[0m    : " << __FILE__ << std::endl; \
    std::cerr << "\033[1mFUNCTION\033[0m: " << __func__ << std::endl; \
    std::cerr << "\033[1mLINE\033[0m    : " << __LINE__ << std::endl; \
    std::cerr << "\033[1mMESSAGE\033[0m : " << m << std::endl; \
    std::cerr << std::endl; \
    exit(2);\
    }

#ifndef NDEBUG
#define ASSERT(c, m) \
    {\
    if (!(c)) { \
        std::cerr << std::endl; \
        std::cerr << "\033[41m\033[1m[ASSERT ERROR]\033[49m\033[0m " << std::endl; \
        std::cerr << "\033[1mFILE\033[0m    : " << __FILE__ << std::endl; \
        std::cerr << "\033[1mFUNCTION\033[0m: " << __func__ << std::endl; \
        std::cerr << "\033[1mLINE\033[0m    : " << __LINE__ << std::endl; \
        std::cerr << "\033[1mMESSAGE\033[0m : " << m << std::endl; \
        std::cerr << std::endl; \
        exit(3); \
    }\
    }
#else
#define ASSERT(c, m) \
    {\
        do \
        { \
        } while(0) \
    }
#endif

#define PRINT_TITLE() \
    std::cout << \
        "fluid-simulation - Tristan Marrec 2021 - Supervised by Ryoichi Ando"\
        << std::endl;
