// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME Dictionary

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
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

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "AliAnalysisTaskDiffCrossSectionsMM.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_AliAnalysisTaskDiffCrossSectionsMM(void *p = 0);
   static void *newArray_AliAnalysisTaskDiffCrossSectionsMM(Long_t size, void *p);
   static void delete_AliAnalysisTaskDiffCrossSectionsMM(void *p);
   static void deleteArray_AliAnalysisTaskDiffCrossSectionsMM(void *p);
   static void destruct_AliAnalysisTaskDiffCrossSectionsMM(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::AliAnalysisTaskDiffCrossSectionsMM*)
   {
      ::AliAnalysisTaskDiffCrossSectionsMM *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::AliAnalysisTaskDiffCrossSectionsMM >(0);
      static ::ROOT::TGenericClassInfo 
         instance("AliAnalysisTaskDiffCrossSectionsMM", ::AliAnalysisTaskDiffCrossSectionsMM::Class_Version(), "AliAnalysisTaskDiffCrossSectionsMM.h", 45,
                  typeid(::AliAnalysisTaskDiffCrossSectionsMM), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::AliAnalysisTaskDiffCrossSectionsMM::Dictionary, isa_proxy, 4,
                  sizeof(::AliAnalysisTaskDiffCrossSectionsMM) );
      instance.SetNew(&new_AliAnalysisTaskDiffCrossSectionsMM);
      instance.SetNewArray(&newArray_AliAnalysisTaskDiffCrossSectionsMM);
      instance.SetDelete(&delete_AliAnalysisTaskDiffCrossSectionsMM);
      instance.SetDeleteArray(&deleteArray_AliAnalysisTaskDiffCrossSectionsMM);
      instance.SetDestructor(&destruct_AliAnalysisTaskDiffCrossSectionsMM);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::AliAnalysisTaskDiffCrossSectionsMM*)
   {
      return GenerateInitInstanceLocal((::AliAnalysisTaskDiffCrossSectionsMM*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::AliAnalysisTaskDiffCrossSectionsMM*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static TClass *AliAnalysisTaskDiffCrossSectionsMMcLcLEventInfo_Dictionary();
   static void AliAnalysisTaskDiffCrossSectionsMMcLcLEventInfo_TClassManip(TClass*);
   static void *new_AliAnalysisTaskDiffCrossSectionsMMcLcLEventInfo(void *p = 0);
   static void *newArray_AliAnalysisTaskDiffCrossSectionsMMcLcLEventInfo(Long_t size, void *p);
   static void delete_AliAnalysisTaskDiffCrossSectionsMMcLcLEventInfo(void *p);
   static void deleteArray_AliAnalysisTaskDiffCrossSectionsMMcLcLEventInfo(void *p);
   static void destruct_AliAnalysisTaskDiffCrossSectionsMMcLcLEventInfo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::AliAnalysisTaskDiffCrossSectionsMM::EventInfo*)
   {
      ::AliAnalysisTaskDiffCrossSectionsMM::EventInfo *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::AliAnalysisTaskDiffCrossSectionsMM::EventInfo));
      static ::ROOT::TGenericClassInfo 
         instance("AliAnalysisTaskDiffCrossSectionsMM::EventInfo", "AliAnalysisTaskDiffCrossSectionsMM.h", 82,
                  typeid(::AliAnalysisTaskDiffCrossSectionsMM::EventInfo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &AliAnalysisTaskDiffCrossSectionsMMcLcLEventInfo_Dictionary, isa_proxy, 4,
                  sizeof(::AliAnalysisTaskDiffCrossSectionsMM::EventInfo) );
      instance.SetNew(&new_AliAnalysisTaskDiffCrossSectionsMMcLcLEventInfo);
      instance.SetNewArray(&newArray_AliAnalysisTaskDiffCrossSectionsMMcLcLEventInfo);
      instance.SetDelete(&delete_AliAnalysisTaskDiffCrossSectionsMMcLcLEventInfo);
      instance.SetDeleteArray(&deleteArray_AliAnalysisTaskDiffCrossSectionsMMcLcLEventInfo);
      instance.SetDestructor(&destruct_AliAnalysisTaskDiffCrossSectionsMMcLcLEventInfo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::AliAnalysisTaskDiffCrossSectionsMM::EventInfo*)
   {
      return GenerateInitInstanceLocal((::AliAnalysisTaskDiffCrossSectionsMM::EventInfo*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::AliAnalysisTaskDiffCrossSectionsMM::EventInfo*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *AliAnalysisTaskDiffCrossSectionsMMcLcLEventInfo_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::AliAnalysisTaskDiffCrossSectionsMM::EventInfo*)0x0)->GetClass();
      AliAnalysisTaskDiffCrossSectionsMMcLcLEventInfo_TClassManip(theClass);
   return theClass;
   }

   static void AliAnalysisTaskDiffCrossSectionsMMcLcLEventInfo_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *AliAnalysisTaskDiffCrossSectionsMMcLcLADV0_Dictionary();
   static void AliAnalysisTaskDiffCrossSectionsMMcLcLADV0_TClassManip(TClass*);
   static void *new_AliAnalysisTaskDiffCrossSectionsMMcLcLADV0(void *p = 0);
   static void *newArray_AliAnalysisTaskDiffCrossSectionsMMcLcLADV0(Long_t size, void *p);
   static void delete_AliAnalysisTaskDiffCrossSectionsMMcLcLADV0(void *p);
   static void deleteArray_AliAnalysisTaskDiffCrossSectionsMMcLcLADV0(void *p);
   static void destruct_AliAnalysisTaskDiffCrossSectionsMMcLcLADV0(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::AliAnalysisTaskDiffCrossSectionsMM::ADV0*)
   {
      ::AliAnalysisTaskDiffCrossSectionsMM::ADV0 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::AliAnalysisTaskDiffCrossSectionsMM::ADV0));
      static ::ROOT::TGenericClassInfo 
         instance("AliAnalysisTaskDiffCrossSectionsMM::ADV0", "AliAnalysisTaskDiffCrossSectionsMM.h", 116,
                  typeid(::AliAnalysisTaskDiffCrossSectionsMM::ADV0), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &AliAnalysisTaskDiffCrossSectionsMMcLcLADV0_Dictionary, isa_proxy, 4,
                  sizeof(::AliAnalysisTaskDiffCrossSectionsMM::ADV0) );
      instance.SetNew(&new_AliAnalysisTaskDiffCrossSectionsMMcLcLADV0);
      instance.SetNewArray(&newArray_AliAnalysisTaskDiffCrossSectionsMMcLcLADV0);
      instance.SetDelete(&delete_AliAnalysisTaskDiffCrossSectionsMMcLcLADV0);
      instance.SetDeleteArray(&deleteArray_AliAnalysisTaskDiffCrossSectionsMMcLcLADV0);
      instance.SetDestructor(&destruct_AliAnalysisTaskDiffCrossSectionsMMcLcLADV0);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::AliAnalysisTaskDiffCrossSectionsMM::ADV0*)
   {
      return GenerateInitInstanceLocal((::AliAnalysisTaskDiffCrossSectionsMM::ADV0*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::AliAnalysisTaskDiffCrossSectionsMM::ADV0*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *AliAnalysisTaskDiffCrossSectionsMMcLcLADV0_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::AliAnalysisTaskDiffCrossSectionsMM::ADV0*)0x0)->GetClass();
      AliAnalysisTaskDiffCrossSectionsMMcLcLADV0_TClassManip(theClass);
   return theClass;
   }

   static void AliAnalysisTaskDiffCrossSectionsMMcLcLADV0_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *AliAnalysisTaskDiffCrossSectionsMMcLcLZDC_Dictionary();
   static void AliAnalysisTaskDiffCrossSectionsMMcLcLZDC_TClassManip(TClass*);
   static void *new_AliAnalysisTaskDiffCrossSectionsMMcLcLZDC(void *p = 0);
   static void *newArray_AliAnalysisTaskDiffCrossSectionsMMcLcLZDC(Long_t size, void *p);
   static void delete_AliAnalysisTaskDiffCrossSectionsMMcLcLZDC(void *p);
   static void deleteArray_AliAnalysisTaskDiffCrossSectionsMMcLcLZDC(void *p);
   static void destruct_AliAnalysisTaskDiffCrossSectionsMMcLcLZDC(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::AliAnalysisTaskDiffCrossSectionsMM::ZDC*)
   {
      ::AliAnalysisTaskDiffCrossSectionsMM::ZDC *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::AliAnalysisTaskDiffCrossSectionsMM::ZDC));
      static ::ROOT::TGenericClassInfo 
         instance("AliAnalysisTaskDiffCrossSectionsMM::ZDC", "AliAnalysisTaskDiffCrossSectionsMM.h", 146,
                  typeid(::AliAnalysisTaskDiffCrossSectionsMM::ZDC), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &AliAnalysisTaskDiffCrossSectionsMMcLcLZDC_Dictionary, isa_proxy, 4,
                  sizeof(::AliAnalysisTaskDiffCrossSectionsMM::ZDC) );
      instance.SetNew(&new_AliAnalysisTaskDiffCrossSectionsMMcLcLZDC);
      instance.SetNewArray(&newArray_AliAnalysisTaskDiffCrossSectionsMMcLcLZDC);
      instance.SetDelete(&delete_AliAnalysisTaskDiffCrossSectionsMMcLcLZDC);
      instance.SetDeleteArray(&deleteArray_AliAnalysisTaskDiffCrossSectionsMMcLcLZDC);
      instance.SetDestructor(&destruct_AliAnalysisTaskDiffCrossSectionsMMcLcLZDC);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::AliAnalysisTaskDiffCrossSectionsMM::ZDC*)
   {
      return GenerateInitInstanceLocal((::AliAnalysisTaskDiffCrossSectionsMM::ZDC*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::AliAnalysisTaskDiffCrossSectionsMM::ZDC*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *AliAnalysisTaskDiffCrossSectionsMMcLcLZDC_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::AliAnalysisTaskDiffCrossSectionsMM::ZDC*)0x0)->GetClass();
      AliAnalysisTaskDiffCrossSectionsMMcLcLZDC_TClassManip(theClass);
   return theClass;
   }

   static void AliAnalysisTaskDiffCrossSectionsMMcLcLZDC_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *AliAnalysisTaskDiffCrossSectionsMMcLcLVtxInfo_Dictionary();
   static void AliAnalysisTaskDiffCrossSectionsMMcLcLVtxInfo_TClassManip(TClass*);
   static void *new_AliAnalysisTaskDiffCrossSectionsMMcLcLVtxInfo(void *p = 0);
   static void *newArray_AliAnalysisTaskDiffCrossSectionsMMcLcLVtxInfo(Long_t size, void *p);
   static void delete_AliAnalysisTaskDiffCrossSectionsMMcLcLVtxInfo(void *p);
   static void deleteArray_AliAnalysisTaskDiffCrossSectionsMMcLcLVtxInfo(void *p);
   static void destruct_AliAnalysisTaskDiffCrossSectionsMMcLcLVtxInfo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::AliAnalysisTaskDiffCrossSectionsMM::VtxInfo*)
   {
      ::AliAnalysisTaskDiffCrossSectionsMM::VtxInfo *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::AliAnalysisTaskDiffCrossSectionsMM::VtxInfo));
      static ::ROOT::TGenericClassInfo 
         instance("AliAnalysisTaskDiffCrossSectionsMM::VtxInfo", "AliAnalysisTaskDiffCrossSectionsMM.h", 173,
                  typeid(::AliAnalysisTaskDiffCrossSectionsMM::VtxInfo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &AliAnalysisTaskDiffCrossSectionsMMcLcLVtxInfo_Dictionary, isa_proxy, 4,
                  sizeof(::AliAnalysisTaskDiffCrossSectionsMM::VtxInfo) );
      instance.SetNew(&new_AliAnalysisTaskDiffCrossSectionsMMcLcLVtxInfo);
      instance.SetNewArray(&newArray_AliAnalysisTaskDiffCrossSectionsMMcLcLVtxInfo);
      instance.SetDelete(&delete_AliAnalysisTaskDiffCrossSectionsMMcLcLVtxInfo);
      instance.SetDeleteArray(&deleteArray_AliAnalysisTaskDiffCrossSectionsMMcLcLVtxInfo);
      instance.SetDestructor(&destruct_AliAnalysisTaskDiffCrossSectionsMMcLcLVtxInfo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::AliAnalysisTaskDiffCrossSectionsMM::VtxInfo*)
   {
      return GenerateInitInstanceLocal((::AliAnalysisTaskDiffCrossSectionsMM::VtxInfo*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::AliAnalysisTaskDiffCrossSectionsMM::VtxInfo*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *AliAnalysisTaskDiffCrossSectionsMMcLcLVtxInfo_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::AliAnalysisTaskDiffCrossSectionsMM::VtxInfo*)0x0)->GetClass();
      AliAnalysisTaskDiffCrossSectionsMMcLcLVtxInfo_TClassManip(theClass);
   return theClass;
   }

   static void AliAnalysisTaskDiffCrossSectionsMMcLcLVtxInfo_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_AliAnalysisTaskDiffCrossSectionsMMcLcLTreeData(void *p = 0);
   static void *newArray_AliAnalysisTaskDiffCrossSectionsMMcLcLTreeData(Long_t size, void *p);
   static void delete_AliAnalysisTaskDiffCrossSectionsMMcLcLTreeData(void *p);
   static void deleteArray_AliAnalysisTaskDiffCrossSectionsMMcLcLTreeData(void *p);
   static void destruct_AliAnalysisTaskDiffCrossSectionsMMcLcLTreeData(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::AliAnalysisTaskDiffCrossSectionsMM::TreeData*)
   {
      ::AliAnalysisTaskDiffCrossSectionsMM::TreeData *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::AliAnalysisTaskDiffCrossSectionsMM::TreeData >(0);
      static ::ROOT::TGenericClassInfo 
         instance("AliAnalysisTaskDiffCrossSectionsMM::TreeData", ::AliAnalysisTaskDiffCrossSectionsMM::TreeData::Class_Version(), "AliAnalysisTaskDiffCrossSectionsMM.h", 187,
                  typeid(::AliAnalysisTaskDiffCrossSectionsMM::TreeData), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::AliAnalysisTaskDiffCrossSectionsMM::TreeData::Dictionary, isa_proxy, 4,
                  sizeof(::AliAnalysisTaskDiffCrossSectionsMM::TreeData) );
      instance.SetNew(&new_AliAnalysisTaskDiffCrossSectionsMMcLcLTreeData);
      instance.SetNewArray(&newArray_AliAnalysisTaskDiffCrossSectionsMMcLcLTreeData);
      instance.SetDelete(&delete_AliAnalysisTaskDiffCrossSectionsMMcLcLTreeData);
      instance.SetDeleteArray(&deleteArray_AliAnalysisTaskDiffCrossSectionsMMcLcLTreeData);
      instance.SetDestructor(&destruct_AliAnalysisTaskDiffCrossSectionsMMcLcLTreeData);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::AliAnalysisTaskDiffCrossSectionsMM::TreeData*)
   {
      return GenerateInitInstanceLocal((::AliAnalysisTaskDiffCrossSectionsMM::TreeData*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::AliAnalysisTaskDiffCrossSectionsMM::TreeData*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfo(void *p = 0);
   static void *newArray_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfo(Long_t size, void *p);
   static void delete_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfo(void *p);
   static void deleteArray_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfo(void *p);
   static void destruct_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfo(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo*)
   {
      ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo >(0);
      static ::ROOT::TGenericClassInfo 
         instance("AliAnalysisTaskDiffCrossSectionsMM::MCInfo", ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::Class_Version(), "AliAnalysisTaskDiffCrossSectionsMM.h", 216,
                  typeid(::AliAnalysisTaskDiffCrossSectionsMM::MCInfo), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::Dictionary, isa_proxy, 4,
                  sizeof(::AliAnalysisTaskDiffCrossSectionsMM::MCInfo) );
      instance.SetNew(&new_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfo);
      instance.SetNewArray(&newArray_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfo);
      instance.SetDelete(&delete_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfo);
      instance.SetDeleteArray(&deleteArray_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfo);
      instance.SetDestructor(&destruct_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfo);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo*)
   {
      return GenerateInitInstanceLocal((::AliAnalysisTaskDiffCrossSectionsMM::MCInfo*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static TClass *AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDetGenLevel_Dictionary();
   static void AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDetGenLevel_TClassManip(TClass*);
   static void *new_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDetGenLevel(void *p = 0);
   static void *newArray_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDetGenLevel(Long_t size, void *p);
   static void delete_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDetGenLevel(void *p);
   static void deleteArray_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDetGenLevel(void *p);
   static void destruct_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDetGenLevel(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DetGenLevel*)
   {
      ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DetGenLevel *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DetGenLevel));
      static ::ROOT::TGenericClassInfo 
         instance("AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DetGenLevel", "AliAnalysisTaskDiffCrossSectionsMM.h", 229,
                  typeid(::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DetGenLevel), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDetGenLevel_Dictionary, isa_proxy, 4,
                  sizeof(::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DetGenLevel) );
      instance.SetNew(&new_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDetGenLevel);
      instance.SetNewArray(&newArray_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDetGenLevel);
      instance.SetDelete(&delete_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDetGenLevel);
      instance.SetDeleteArray(&deleteArray_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDetGenLevel);
      instance.SetDestructor(&destruct_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDetGenLevel);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DetGenLevel*)
   {
      return GenerateInitInstanceLocal((::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DetGenLevel*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DetGenLevel*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDetGenLevel_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DetGenLevel*)0x0)->GetClass();
      AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDetGenLevel_TClassManip(theClass);
   return theClass;
   }

   static void AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDetGenLevel_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDiffSystem_Dictionary();
   static void AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDiffSystem_TClassManip(TClass*);
   static void *new_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDiffSystem(void *p = 0);
   static void *newArray_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDiffSystem(Long_t size, void *p);
   static void delete_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDiffSystem(void *p);
   static void deleteArray_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDiffSystem(void *p);
   static void destruct_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDiffSystem(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DiffSystem*)
   {
      ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DiffSystem *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DiffSystem));
      static ::ROOT::TGenericClassInfo 
         instance("AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DiffSystem", "AliAnalysisTaskDiffCrossSectionsMM.h", 260,
                  typeid(::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DiffSystem), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDiffSystem_Dictionary, isa_proxy, 4,
                  sizeof(::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DiffSystem) );
      instance.SetNew(&new_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDiffSystem);
      instance.SetNewArray(&newArray_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDiffSystem);
      instance.SetDelete(&delete_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDiffSystem);
      instance.SetDeleteArray(&deleteArray_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDiffSystem);
      instance.SetDestructor(&destruct_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDiffSystem);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DiffSystem*)
   {
      return GenerateInitInstanceLocal((::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DiffSystem*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DiffSystem*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDiffSystem_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DiffSystem*)0x0)->GetClass();
      AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDiffSystem_TClassManip(theClass);
   return theClass;
   }

   static void AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDiffSystem_TClassManip(TClass* ){
   }

} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr AliAnalysisTaskDiffCrossSectionsMM::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *AliAnalysisTaskDiffCrossSectionsMM::Class_Name()
{
   return "AliAnalysisTaskDiffCrossSectionsMM";
}

//______________________________________________________________________________
const char *AliAnalysisTaskDiffCrossSectionsMM::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::AliAnalysisTaskDiffCrossSectionsMM*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int AliAnalysisTaskDiffCrossSectionsMM::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::AliAnalysisTaskDiffCrossSectionsMM*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *AliAnalysisTaskDiffCrossSectionsMM::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::AliAnalysisTaskDiffCrossSectionsMM*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *AliAnalysisTaskDiffCrossSectionsMM::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::AliAnalysisTaskDiffCrossSectionsMM*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr AliAnalysisTaskDiffCrossSectionsMM::TreeData::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *AliAnalysisTaskDiffCrossSectionsMM::TreeData::Class_Name()
{
   return "AliAnalysisTaskDiffCrossSectionsMM::TreeData";
}

//______________________________________________________________________________
const char *AliAnalysisTaskDiffCrossSectionsMM::TreeData::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::AliAnalysisTaskDiffCrossSectionsMM::TreeData*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int AliAnalysisTaskDiffCrossSectionsMM::TreeData::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::AliAnalysisTaskDiffCrossSectionsMM::TreeData*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *AliAnalysisTaskDiffCrossSectionsMM::TreeData::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::AliAnalysisTaskDiffCrossSectionsMM::TreeData*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *AliAnalysisTaskDiffCrossSectionsMM::TreeData::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::AliAnalysisTaskDiffCrossSectionsMM::TreeData*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr AliAnalysisTaskDiffCrossSectionsMM::MCInfo::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *AliAnalysisTaskDiffCrossSectionsMM::MCInfo::Class_Name()
{
   return "AliAnalysisTaskDiffCrossSectionsMM::MCInfo";
}

//______________________________________________________________________________
const char *AliAnalysisTaskDiffCrossSectionsMM::MCInfo::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int AliAnalysisTaskDiffCrossSectionsMM::MCInfo::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *AliAnalysisTaskDiffCrossSectionsMM::MCInfo::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *AliAnalysisTaskDiffCrossSectionsMM::MCInfo::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD2(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void AliAnalysisTaskDiffCrossSectionsMM::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliAnalysisTaskDiffCrossSectionsMM.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(AliAnalysisTaskDiffCrossSectionsMM::Class(),this);
   } else {
      R__b.WriteClassBuffer(AliAnalysisTaskDiffCrossSectionsMM::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_AliAnalysisTaskDiffCrossSectionsMM(void *p) {
      return  p ? new(p) ::AliAnalysisTaskDiffCrossSectionsMM : new ::AliAnalysisTaskDiffCrossSectionsMM;
   }
   static void *newArray_AliAnalysisTaskDiffCrossSectionsMM(Long_t nElements, void *p) {
      return p ? new(p) ::AliAnalysisTaskDiffCrossSectionsMM[nElements] : new ::AliAnalysisTaskDiffCrossSectionsMM[nElements];
   }
   // Wrapper around operator delete
   static void delete_AliAnalysisTaskDiffCrossSectionsMM(void *p) {
      delete ((::AliAnalysisTaskDiffCrossSectionsMM*)p);
   }
   static void deleteArray_AliAnalysisTaskDiffCrossSectionsMM(void *p) {
      delete [] ((::AliAnalysisTaskDiffCrossSectionsMM*)p);
   }
   static void destruct_AliAnalysisTaskDiffCrossSectionsMM(void *p) {
      typedef ::AliAnalysisTaskDiffCrossSectionsMM current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::AliAnalysisTaskDiffCrossSectionsMM

namespace ROOT {
   // Wrappers around operator new
   static void *new_AliAnalysisTaskDiffCrossSectionsMMcLcLEventInfo(void *p) {
      return  p ? new(p) ::AliAnalysisTaskDiffCrossSectionsMM::EventInfo : new ::AliAnalysisTaskDiffCrossSectionsMM::EventInfo;
   }
   static void *newArray_AliAnalysisTaskDiffCrossSectionsMMcLcLEventInfo(Long_t nElements, void *p) {
      return p ? new(p) ::AliAnalysisTaskDiffCrossSectionsMM::EventInfo[nElements] : new ::AliAnalysisTaskDiffCrossSectionsMM::EventInfo[nElements];
   }
   // Wrapper around operator delete
   static void delete_AliAnalysisTaskDiffCrossSectionsMMcLcLEventInfo(void *p) {
      delete ((::AliAnalysisTaskDiffCrossSectionsMM::EventInfo*)p);
   }
   static void deleteArray_AliAnalysisTaskDiffCrossSectionsMMcLcLEventInfo(void *p) {
      delete [] ((::AliAnalysisTaskDiffCrossSectionsMM::EventInfo*)p);
   }
   static void destruct_AliAnalysisTaskDiffCrossSectionsMMcLcLEventInfo(void *p) {
      typedef ::AliAnalysisTaskDiffCrossSectionsMM::EventInfo current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::AliAnalysisTaskDiffCrossSectionsMM::EventInfo

namespace ROOT {
   // Wrappers around operator new
   static void *new_AliAnalysisTaskDiffCrossSectionsMMcLcLADV0(void *p) {
      return  p ? new(p) ::AliAnalysisTaskDiffCrossSectionsMM::ADV0 : new ::AliAnalysisTaskDiffCrossSectionsMM::ADV0;
   }
   static void *newArray_AliAnalysisTaskDiffCrossSectionsMMcLcLADV0(Long_t nElements, void *p) {
      return p ? new(p) ::AliAnalysisTaskDiffCrossSectionsMM::ADV0[nElements] : new ::AliAnalysisTaskDiffCrossSectionsMM::ADV0[nElements];
   }
   // Wrapper around operator delete
   static void delete_AliAnalysisTaskDiffCrossSectionsMMcLcLADV0(void *p) {
      delete ((::AliAnalysisTaskDiffCrossSectionsMM::ADV0*)p);
   }
   static void deleteArray_AliAnalysisTaskDiffCrossSectionsMMcLcLADV0(void *p) {
      delete [] ((::AliAnalysisTaskDiffCrossSectionsMM::ADV0*)p);
   }
   static void destruct_AliAnalysisTaskDiffCrossSectionsMMcLcLADV0(void *p) {
      typedef ::AliAnalysisTaskDiffCrossSectionsMM::ADV0 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::AliAnalysisTaskDiffCrossSectionsMM::ADV0

namespace ROOT {
   // Wrappers around operator new
   static void *new_AliAnalysisTaskDiffCrossSectionsMMcLcLZDC(void *p) {
      return  p ? new(p) ::AliAnalysisTaskDiffCrossSectionsMM::ZDC : new ::AliAnalysisTaskDiffCrossSectionsMM::ZDC;
   }
   static void *newArray_AliAnalysisTaskDiffCrossSectionsMMcLcLZDC(Long_t nElements, void *p) {
      return p ? new(p) ::AliAnalysisTaskDiffCrossSectionsMM::ZDC[nElements] : new ::AliAnalysisTaskDiffCrossSectionsMM::ZDC[nElements];
   }
   // Wrapper around operator delete
   static void delete_AliAnalysisTaskDiffCrossSectionsMMcLcLZDC(void *p) {
      delete ((::AliAnalysisTaskDiffCrossSectionsMM::ZDC*)p);
   }
   static void deleteArray_AliAnalysisTaskDiffCrossSectionsMMcLcLZDC(void *p) {
      delete [] ((::AliAnalysisTaskDiffCrossSectionsMM::ZDC*)p);
   }
   static void destruct_AliAnalysisTaskDiffCrossSectionsMMcLcLZDC(void *p) {
      typedef ::AliAnalysisTaskDiffCrossSectionsMM::ZDC current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::AliAnalysisTaskDiffCrossSectionsMM::ZDC

namespace ROOT {
   // Wrappers around operator new
   static void *new_AliAnalysisTaskDiffCrossSectionsMMcLcLVtxInfo(void *p) {
      return  p ? new(p) ::AliAnalysisTaskDiffCrossSectionsMM::VtxInfo : new ::AliAnalysisTaskDiffCrossSectionsMM::VtxInfo;
   }
   static void *newArray_AliAnalysisTaskDiffCrossSectionsMMcLcLVtxInfo(Long_t nElements, void *p) {
      return p ? new(p) ::AliAnalysisTaskDiffCrossSectionsMM::VtxInfo[nElements] : new ::AliAnalysisTaskDiffCrossSectionsMM::VtxInfo[nElements];
   }
   // Wrapper around operator delete
   static void delete_AliAnalysisTaskDiffCrossSectionsMMcLcLVtxInfo(void *p) {
      delete ((::AliAnalysisTaskDiffCrossSectionsMM::VtxInfo*)p);
   }
   static void deleteArray_AliAnalysisTaskDiffCrossSectionsMMcLcLVtxInfo(void *p) {
      delete [] ((::AliAnalysisTaskDiffCrossSectionsMM::VtxInfo*)p);
   }
   static void destruct_AliAnalysisTaskDiffCrossSectionsMMcLcLVtxInfo(void *p) {
      typedef ::AliAnalysisTaskDiffCrossSectionsMM::VtxInfo current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::AliAnalysisTaskDiffCrossSectionsMM::VtxInfo

//______________________________________________________________________________
void AliAnalysisTaskDiffCrossSectionsMM::TreeData::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliAnalysisTaskDiffCrossSectionsMM::TreeData.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(AliAnalysisTaskDiffCrossSectionsMM::TreeData::Class(),this);
   } else {
      R__b.WriteClassBuffer(AliAnalysisTaskDiffCrossSectionsMM::TreeData::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_AliAnalysisTaskDiffCrossSectionsMMcLcLTreeData(void *p) {
      return  p ? new(p) ::AliAnalysisTaskDiffCrossSectionsMM::TreeData : new ::AliAnalysisTaskDiffCrossSectionsMM::TreeData;
   }
   static void *newArray_AliAnalysisTaskDiffCrossSectionsMMcLcLTreeData(Long_t nElements, void *p) {
      return p ? new(p) ::AliAnalysisTaskDiffCrossSectionsMM::TreeData[nElements] : new ::AliAnalysisTaskDiffCrossSectionsMM::TreeData[nElements];
   }
   // Wrapper around operator delete
   static void delete_AliAnalysisTaskDiffCrossSectionsMMcLcLTreeData(void *p) {
      delete ((::AliAnalysisTaskDiffCrossSectionsMM::TreeData*)p);
   }
   static void deleteArray_AliAnalysisTaskDiffCrossSectionsMMcLcLTreeData(void *p) {
      delete [] ((::AliAnalysisTaskDiffCrossSectionsMM::TreeData*)p);
   }
   static void destruct_AliAnalysisTaskDiffCrossSectionsMMcLcLTreeData(void *p) {
      typedef ::AliAnalysisTaskDiffCrossSectionsMM::TreeData current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::AliAnalysisTaskDiffCrossSectionsMM::TreeData

//______________________________________________________________________________
void AliAnalysisTaskDiffCrossSectionsMM::MCInfo::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliAnalysisTaskDiffCrossSectionsMM::MCInfo.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(AliAnalysisTaskDiffCrossSectionsMM::MCInfo::Class(),this);
   } else {
      R__b.WriteClassBuffer(AliAnalysisTaskDiffCrossSectionsMM::MCInfo::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfo(void *p) {
      return  p ? new(p) ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo : new ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo;
   }
   static void *newArray_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfo(Long_t nElements, void *p) {
      return p ? new(p) ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo[nElements] : new ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo[nElements];
   }
   // Wrapper around operator delete
   static void delete_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfo(void *p) {
      delete ((::AliAnalysisTaskDiffCrossSectionsMM::MCInfo*)p);
   }
   static void deleteArray_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfo(void *p) {
      delete [] ((::AliAnalysisTaskDiffCrossSectionsMM::MCInfo*)p);
   }
   static void destruct_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfo(void *p) {
      typedef ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo

namespace ROOT {
   // Wrappers around operator new
   static void *new_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDetGenLevel(void *p) {
      return  p ? new(p) ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DetGenLevel : new ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DetGenLevel;
   }
   static void *newArray_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDetGenLevel(Long_t nElements, void *p) {
      return p ? new(p) ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DetGenLevel[nElements] : new ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DetGenLevel[nElements];
   }
   // Wrapper around operator delete
   static void delete_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDetGenLevel(void *p) {
      delete ((::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DetGenLevel*)p);
   }
   static void deleteArray_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDetGenLevel(void *p) {
      delete [] ((::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DetGenLevel*)p);
   }
   static void destruct_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDetGenLevel(void *p) {
      typedef ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DetGenLevel current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DetGenLevel

namespace ROOT {
   // Wrappers around operator new
   static void *new_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDiffSystem(void *p) {
      return  p ? new(p) ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DiffSystem : new ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DiffSystem;
   }
   static void *newArray_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDiffSystem(Long_t nElements, void *p) {
      return p ? new(p) ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DiffSystem[nElements] : new ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DiffSystem[nElements];
   }
   // Wrapper around operator delete
   static void delete_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDiffSystem(void *p) {
      delete ((::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DiffSystem*)p);
   }
   static void deleteArray_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDiffSystem(void *p) {
      delete [] ((::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DiffSystem*)p);
   }
   static void destruct_AliAnalysisTaskDiffCrossSectionsMMcLcLMCInfocLcLDiffSystem(void *p) {
      typedef ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DiffSystem current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DiffSystem

namespace {
  void TriggerDictionaryInitialization_Dictionary_Impl() {
    static const char* headers[] = {
"AliAnalysisTaskDiffCrossSectionsMM.h",
0
    };
    static const char* includePaths[] = {
"/home/user/alice/sw/ubuntu1804_x86-64/AliRoot/latest/include",
"/home/user/alice/sw/ubuntu1804_x86-64/AliPhysics/latest/include",
"/home/user/alice/sw/ubuntu1804_x86-64/ROOT/v6-10-08-1/include",
"/home/user/cernbox/ALICE/Diffractive-Combinatorics/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "Dictionary dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate("$clingAutoload$AliAnalysisTaskDiffCrossSectionsMM.h")))  AliAnalysisTaskDiffCrossSectionsMM;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "Dictionary dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "AliAnalysisTaskDiffCrossSectionsMM.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"AliAnalysisTaskDiffCrossSectionsMM", payloadCode, "@",
"AliAnalysisTaskDiffCrossSectionsMM::ADV0", payloadCode, "@",
"AliAnalysisTaskDiffCrossSectionsMM::EventInfo", payloadCode, "@",
"AliAnalysisTaskDiffCrossSectionsMM::MCInfo", payloadCode, "@",
"AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DetGenLevel", payloadCode, "@",
"AliAnalysisTaskDiffCrossSectionsMM::MCInfo::DiffSystem", payloadCode, "@",
"AliAnalysisTaskDiffCrossSectionsMM::TreeData", payloadCode, "@",
"AliAnalysisTaskDiffCrossSectionsMM::VtxInfo", payloadCode, "@",
"AliAnalysisTaskDiffCrossSectionsMM::ZDC", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("Dictionary",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_Dictionary_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_Dictionary_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_Dictionary() {
  TriggerDictionaryInitialization_Dictionary_Impl();
}
