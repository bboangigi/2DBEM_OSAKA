SET(ITEM_PREFIX "LIB")
SET(ITEM_NAME   "Libraries_2DBEM_OSAKA")

## SRC 폴더 내의 Library 호출 목록임. 순서는 중요하지 않음
## (1) 만약 library 호출 폴더 이름이 바뀌게 된다면?  폴더명(윈도우탐색기): libe2DBEM -> lib2DBEM_OSAKA
##  1) src 폴더 내의 CMakeLists.txt: lib2DBEM -> lib2DBEM_01 각 library 연결 폴더 이름 전부 바꾸기! 
##  2) 상위 폴더 내의 CMakeLists.txt: Library call folder 란에 set(LIB_NAME_lib2DBEM lib2DBEM) -> set(LIB_NAME_lib2DBEM lib2DBEM_01)
##  3) 상위 폴더 내의 CMakeLists.txt: Application 란에 set(APP_DEPEND_app2DBEM lib2DBEM) -> set(APP_DEPEND_app2DBEM lib2DBEM_01)
##  4) src 폴더 내의 CMakeLists.txt: ITEM_PREFIX 및 ITEM_NAME은 절대로 바꾸지 말 것!

## SRC 폴더 내의 ITEM_NAME 바꾸는 절차임.
## (2) 만약 ITEM_NAME이 lib2DBEM에서 Libraries_2DBEM_OSAK로 바귀게 된다면?
##  1) scr 폴더 내의 CMakeLists.txt: lib2DBEM -> Libraries_2DBEM_OSAKA 각 library 명칭 (폴더명이랑 상관 없음) 바꾸기!
##  2) 상위 폴더 내의 CMakeLists.txt: set(LIB_NAME_lib2DBEM lib2DBEM_OSAKA) -> set(LIB_NAME_Libraries_2DBEM_OSAKA lib2DBEM_OSAKA)

set(${ITEM_PREFIX}_SRC_${ITEM_NAME}
    lib2DBEM_OSAKA/glb2DBEM.F90
    lib2DBEM_OSAKA/pkg2DBEM_OSAKA.F90

    lib2DBEM_OSAKA/MOTION_M.F90
    lib2DBEM_OSAKA/FORCE_M.F90
    lib2DBEM_OSAKA/OFFSET_M.F90
    lib2DBEM_OSAKA/SOLVE_M.F90
)

##... Compile
AutoCompile(
    IS_COMPILE_LIB      True
    IS_DYNAMIC_LIB      "${COMPILE_DYNAMIC_LIB}"
    INSTALL_NAME        "${${ITEM_PREFIX}_NAME_${ITEM_NAME}}"
    INSTALL_SRC         "${${ITEM_PREFIX}_SRC_${ITEM_NAME}}"
    INSTALL_PATH        "${${ITEM_PREFIX}_PATH_${ITEM_NAME}}"
    INSTALL_DEPEND      "${${ITEM_PREFIX}_DEPEND_${ITEM_NAME}}"
    INSTALL_EXTLIB      "${${ITEM_PREFIX}_EXTLIB_${ITEM_NAME}}"
    INSTALL_EXTINC      "${${ITEM_PREFIX}_EXTINC_${ITEM_NAME}}"
    INSTALL_EXTFLAG     "${${ITEM_PREFIX}_EXTFLAG_${ITEM_NAME}}"
)