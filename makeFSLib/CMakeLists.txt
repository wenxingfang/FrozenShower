
set(LCRelation_check_src ./src/LCRelation_check/main.cpp ./src/LCRelation_check/SLCIORdr.cpp ./src/utils/my_utils.cpp)
set(LCRelation_check_inc ./src/LCRelation_check/ ./src/utils/)

set(main_src ./src/FS/main.cpp ./src/FS/SLCIORdr.cpp ./src/utils/my_utils.cpp)
set(main_inc ./src/FS/ ./src/utils/)

set(src_FS_digi ./src/FS_digi/main.cpp ./src/FS_digi/SLCIORdr.cpp ./src/utils/my_utils.cpp)
set(inc_FS_digi ./src/FS_digi/ ./src/utils/)

set(src_Ecal_check ./src/ECAL_check/main.cpp)
#set(main_src_1 ./src/makeLib/main.cpp ./src/makeLib/makeLib.cpp ./src/utils/my_utils.cpp)
#set(main_inc_1 ./src/makeLib/ ./src/utils/)
#set(main_src_1 ./src/makeLib_v1/main.cpp ./src/makeLib_v1/makeLib.cpp ./src/utils/my_utils.cpp)
#set(main_inc_1 ./src/makeLib_v1/ ./src/utils/)
#set(main_src_1 ./src/makeLib_v2/main.cpp ./src/makeLib_v2/makeLib.cpp ./src/utils/my_utils.cpp)
#set(main_inc_1 ./src/makeLib_v2/ ./src/utils/)
#set(main_src_1 ./src/makeLib_v3/main.cpp ./src/makeLib_v3/makeLib.cpp ./src/utils/my_utils.cpp)
#set(main_inc_1 ./src/makeLib_v3/ ./src/utils/)
#set(main_src_1 ./src/makeLib_v4/main.cpp ./src/makeLib_v4/makeLib.cpp ./src/utils/my_utils.cpp)
#set(main_inc_1 ./src/makeLib_v4/ ./src/utils/)
#set(main_src_1 ./src/makeLib_v5/main.cpp ./src/makeLib_v5/makeLib.cpp ./src/utils/my_utils.cpp)
#set(main_inc_1 ./src/makeLib_v5/ ./src/utils/)
#set(main_src_1 ./src/makeLib_v6/main.cpp ./src/makeLib_v6/makeLib.cpp ./src/utils/my_utils.cpp)
#set(main_inc_1 ./src/makeLib_v6/ ./src/utils/)
#set(main_src_1 ./src/makeLib_sp/main.cpp ./src/makeLib_sp/makeLib.cpp ./src/utils/my_utils.cpp)
#set(main_inc_1 ./src/makeLib_sp/ ./src/utils/)
set(main_src_1 ./src/makeLib_sp_v1/main.cpp ./src/makeLib_sp_v1/makeLib.cpp ./src/utils/my_utils.cpp)
set(main_inc_1 ./src/makeLib_sp_v1/ ./src/utils/)

set(main_src_2 ./src/checkLib_sp/main.cpp ./src/checkLib_sp/checkLib.cpp ./src/utils/my_utils.cpp)
set(main_inc_2 ./src/checkLib_sp/ ./src/utils/)

set(main_src_3 ./src/makeLib_sp_faiss/main.cpp ./src/makeLib_sp_faiss/makeLib.cpp ./src/utils/my_utils.cpp)
set(main_inc_3 ./src/makeLib_sp_faiss/ ./src/utils/)

set(main_src_4 ./src/writeIndex_faiss/main.cpp)
set(main_src_5 ./src/readIndex_faiss/main.cpp)

set(src_txt ./src/makeTxtRoot/main.cpp ./src/makeTxtRoot/makeLib.cpp ./src/utils/my_utils.cpp)
set(inc_txt ./src/makeTxtRoot/ ./src/utils/)


SET(CMAKE_MODULE_PATH /cvmfs/cepc.ihep.ac.cn/software/cepcsoft/x86_64-sl6-gcc49/external/ROOT/5.34.18/install/etc/cmake)

find_package(LCIO)
#find_package(ROOT COMPONENTS RIO Tree TreePlayer MathCore Net Graf3d Graf Gpad REQUIRED)
find_package(ROOT)

if(ROOT_FOUND)
    message("found ROOT: ${ROOT_INCLUDE_DIR} ${ROOT_LIBRARIES}")
endif(ROOT_FOUND)

if(LCIO_FOUND)
    message("found LCIO: ${LCIO_INCLUDE_DIRS} ${LCIO_LIBRARIES}")
endif(LCIO_FOUND)

INCLUDE_DIRECTORIES(${LCRelation_check_inc})
INCLUDE_DIRECTORIES(${inc_FS_digi})
INCLUDE_DIRECTORIES(${inc_txt})
INCLUDE_DIRECTORIES(${main_inc})
INCLUDE_DIRECTORIES(${main_inc_1})
INCLUDE_DIRECTORIES(${main_inc_2})
INCLUDE_DIRECTORIES(${main_inc_3})
INCLUDE_DIRECTORIES(${LCIO_INCLUDE_DIRS})
INCLUDE_DIRECTORIES(${ROOT_INCLUDE_DIR})

#############################################################################################################################
set(faiss_inc /junofs/users/wxfang/faiss)
INCLUDE_DIRECTORIES(${faiss_inc})
LINK_LIBRARIES("/junofs/users/wxfang/faiss/libfaiss.so")
ADD_DEFINITIONS( "-std=c++11" )##wxfang
#############################################################################################################################

#ADD_EXECUTABLE(LCRelation_check ${LCRelation_check_src})
#TARGET_LINK_LIBRARIES(LCRelation_check ${LCIO_LIBRARIES} ${ROOT_LIBRARIES})

#ADD_EXECUTABLE(processFS ${main_src})
#TARGET_LINK_LIBRARIES(processFS ${LCIO_LIBRARIES} ${ROOT_LIBRARIES})

ADD_EXECUTABLE(processFS_digi ${src_FS_digi})
TARGET_LINK_LIBRARIES(processFS_digi ${LCIO_LIBRARIES} ${ROOT_LIBRARIES})

#ADD_EXECUTABLE(Ecal_check ${src_Ecal_check})
#TARGET_LINK_LIBRARIES(Ecal_check ${LCIO_LIBRARIES} ${ROOT_LIBRARIES})

#ADD_EXECUTABLE(processMakeSP ${main_src_1})
#TARGET_LINK_LIBRARIES(processMakeSP ${LCIO_LIBRARIES} ${ROOT_LIBRARIES} /junofs/users/wxfang/CEPC/CEPCOFF/FS_study/AutoDict_std__vector__std__vector_int____cxx.so /junofs/users/wxfang/CEPC/CEPCOFF/FS_study/AutoDict_std__vector__std__vector_float____cxx.so)


#ADD_EXECUTABLE(processCheckLibSP ${main_src_2})
#TARGET_LINK_LIBRARIES(processCheckLibSP ${LCIO_LIBRARIES} ${ROOT_LIBRARIES} /junofs/users/wxfang/CEPC/CEPCOFF/FS_study/AutoDict_std__vector__std__vector_int____cxx.so /junofs/users/wxfang/CEPC/CEPCOFF/FS_study/AutoDict_std__vector__std__vector_float____cxx.so)

#ADD_EXECUTABLE(processMakeFaiss ${main_src_3})
#TARGET_LINK_LIBRARIES(processMakeFaiss ${LCIO_LIBRARIES} ${ROOT_LIBRARIES} /junofs/users/wxfang/CEPC/CEPCOFF/FS_study/AutoDict_std__vector__std__vector_int____cxx.so /junofs/users/wxfang/CEPC/CEPCOFF/FS_study/AutoDict_std__vector__std__vector_float____cxx.so)

#ADD_EXECUTABLE(processWriteIndexFaiss ${main_src_4})
#ADD_EXECUTABLE(processReadIndexFaiss ${main_src_5})

#ADD_EXECUTABLE(processTxt ${src_txt})
#TARGET_LINK_LIBRARIES(processTxt ${LCIO_LIBRARIES} ${ROOT_LIBRARIES} /junofs/users/wxfang/CEPC/CEPCOFF/FS_study/AutoDict_std__vector__std__vector_int____cxx.so /junofs/users/wxfang/CEPC/CEPCOFF/FS_study/AutoDict_std__vector__std__vector_float____cxx.so)

