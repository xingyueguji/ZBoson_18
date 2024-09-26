// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME myDict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *vectorlEvectorlEunsignedsPchargRsPgR_Dictionary();
   static void vectorlEvectorlEunsignedsPchargRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEunsignedsPchargRsPgR(void *p = nullptr);
   static void *newArray_vectorlEvectorlEunsignedsPchargRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEunsignedsPchargRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEunsignedsPchargRsPgR(void *p);
   static void destruct_vectorlEvectorlEunsignedsPchargRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<unsigned char> >*)
   {
      vector<vector<unsigned char> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<unsigned char> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<unsigned char> >", -2, "vector", 389,
                  typeid(vector<vector<unsigned char> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEunsignedsPchargRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<unsigned char> >) );
      instance.SetNew(&new_vectorlEvectorlEunsignedsPchargRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEunsignedsPchargRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEunsignedsPchargRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEunsignedsPchargRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEunsignedsPchargRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<unsigned char> > >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<vector<unsigned char> >","std::vector<std::vector<unsigned char, std::allocator<unsigned char> >, std::allocator<std::vector<unsigned char, std::allocator<unsigned char> > > >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<vector<unsigned char> >*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEunsignedsPchargRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<vector<unsigned char> >*>(nullptr))->GetClass();
      vectorlEvectorlEunsignedsPchargRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEunsignedsPchargRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEunsignedsPchargRsPgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<unsigned char> > : new vector<vector<unsigned char> >;
   }
   static void *newArray_vectorlEvectorlEunsignedsPchargRsPgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<vector<unsigned char> >[nElements] : new vector<vector<unsigned char> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEunsignedsPchargRsPgR(void *p) {
      delete (static_cast<vector<vector<unsigned char> >*>(p));
   }
   static void deleteArray_vectorlEvectorlEunsignedsPchargRsPgR(void *p) {
      delete [] (static_cast<vector<vector<unsigned char> >*>(p));
   }
   static void destruct_vectorlEvectorlEunsignedsPchargRsPgR(void *p) {
      typedef vector<vector<unsigned char> > current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<vector<unsigned char> >

namespace ROOT {
   static TClass *vectorlEunsignedsPchargR_Dictionary();
   static void vectorlEunsignedsPchargR_TClassManip(TClass*);
   static void *new_vectorlEunsignedsPchargR(void *p = nullptr);
   static void *newArray_vectorlEunsignedsPchargR(Long_t size, void *p);
   static void delete_vectorlEunsignedsPchargR(void *p);
   static void deleteArray_vectorlEunsignedsPchargR(void *p);
   static void destruct_vectorlEunsignedsPchargR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<unsigned char>*)
   {
      vector<unsigned char> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<unsigned char>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<unsigned char>", -2, "vector", 389,
                  typeid(vector<unsigned char>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEunsignedsPchargR_Dictionary, isa_proxy, 4,
                  sizeof(vector<unsigned char>) );
      instance.SetNew(&new_vectorlEunsignedsPchargR);
      instance.SetNewArray(&newArray_vectorlEunsignedsPchargR);
      instance.SetDelete(&delete_vectorlEunsignedsPchargR);
      instance.SetDeleteArray(&deleteArray_vectorlEunsignedsPchargR);
      instance.SetDestructor(&destruct_vectorlEunsignedsPchargR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<unsigned char> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<unsigned char>","std::vector<unsigned char, std::allocator<unsigned char> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<unsigned char>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEunsignedsPchargR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<unsigned char>*>(nullptr))->GetClass();
      vectorlEunsignedsPchargR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEunsignedsPchargR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEunsignedsPchargR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<unsigned char> : new vector<unsigned char>;
   }
   static void *newArray_vectorlEunsignedsPchargR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<unsigned char>[nElements] : new vector<unsigned char>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEunsignedsPchargR(void *p) {
      delete (static_cast<vector<unsigned char>*>(p));
   }
   static void deleteArray_vectorlEunsignedsPchargR(void *p) {
      delete [] (static_cast<vector<unsigned char>*>(p));
   }
   static void destruct_vectorlEunsignedsPchargR(void *p) {
      typedef vector<unsigned char> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<unsigned char>

namespace {
  void TriggerDictionaryInitialization_myDict_Impl() {
    static const char* headers[] = {
nullptr
    };
    static const char* includePaths[] = {
"/usr/include/root",
"/afs/cern.ch/user/z/zheng/ZBoson_18/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "myDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "myDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("myDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_myDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_myDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_myDict() {
  TriggerDictionaryInitialization_myDict_Impl();
}
