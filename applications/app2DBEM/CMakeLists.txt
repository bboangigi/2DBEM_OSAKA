SET(ITEM_PREFIX "APP")
SET(ITEM_NAME   "app2DBEM")

set(${ITEM_PREFIX}_SRC_${ITEM_NAME}

    src/app2DBEM_OSAKA.F90
)

##... Compile!
## 이쪽은 너무 복잡해서 그대로 default 처리 해서 사용할 것!
AutoCompile(
    IS_COMPILE_EXE      True
    INSTALL_NAME        "${${ITEM_PREFIX}_NAME_${ITEM_NAME}}"
    INSTALL_SRC         "${${ITEM_PREFIX}_SRC_${ITEM_NAME}}"
    INSTALL_PATH        "${${ITEM_PREFIX}_PATH_${ITEM_NAME}}"
    INSTALL_DEPEND      "${${ITEM_PREFIX}_DEPEND_${ITEM_NAME}}"
    INSTALL_EXTLIB      "${${ITEM_PREFIX}_EXTLIB_${ITEM_NAME}}"
    INSTALL_EXTINC      "${${ITEM_PREFIX}_EXTINC_${ITEM_NAME}}"
    INSTALL_EXTFLAG     "${${ITEM_PREFIX}_EXTFLAG_${ITEM_NAME}}"
)
