TEMPLATE = lib

CONFIG -= app_bundle
CONFIG -= qt

macx: {
    INCLUDEPATH += /opt/local/include/eigen3
    INCLUDEPATH += /opt/local/include

    QMAKE_CXXFLAGS += -stdlib=libc++
    LIBS += -stdlib=libc++
} else:unix: {
    INCLUDEPATH += /usr/include/eigen3
}

INCLUDEPATH += include

HEADERS += \
    include/opengnc/estimation/unscented_transform.hpp \
    include/opengnc/estimation/time_update.hpp \
    include/opengnc/estimation/models/measurement/gps.hpp \
    include/opengnc/common/transforms/wgs84.hpp \
    include/opengnc/common/first_order_density.hpp \
    include/opengnc/estimation/measurement_update.hpp \
    include/opengnc/estimation/models/process/rigid_body/constant_acceleration.hpp \
    include/opengnc/estimation/models/process/rigid_body/dwna_covariance_policy.hpp



