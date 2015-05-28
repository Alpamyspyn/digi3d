/****************************************************************************
** Meta object code from reading C++ file 'fs_tcp_client.h'
**
** Created by: The Qt Meta Object Compiler version 67 (Qt 5.4.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../../Exercise6-Package/06-RigidRegistration/fs_tcp_client.h"
#include <QtCore/qbytearray.h>
#include <QtCore/qmetatype.h>
#include <QtCore/QSharedPointer>
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'fs_tcp_client.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 67
#error "This file was generated using the moc from 5.4.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
struct qt_meta_stringdata_FaceShiftClient_t {
    QByteArrayData data[6];
    char stringdata[117];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_FaceShiftClient_t, stringdata) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_FaceShiftClient_t qt_meta_stringdata_FaceShiftClient = {
    {
QT_MOC_LITERAL(0, 0, 15), // "FaceShiftClient"
QT_MOC_LITERAL(1, 16, 16), // "askForConnection"
QT_MOC_LITERAL(2, 33, 0), // ""
QT_MOC_LITERAL(3, 34, 20), // "initializeConnection"
QT_MOC_LITERAL(4, 55, 22), // "receivedNetworkUpdates"
QT_MOC_LITERAL(5, 78, 38) // "QSharedPointer<fs::fsMsgTrack..."

    },
    "FaceShiftClient\0askForConnection\0\0"
    "initializeConnection\0receivedNetworkUpdates\0"
    "QSharedPointer<fs::fsMsgTrackingState>"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_FaceShiftClient[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       3,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    2,   29,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       3,    0,   34,    2, 0x0a /* Public */,
       4,    1,   35,    2, 0x0a /* Public */,

 // signals: parameters
    QMetaType::Void, QMetaType::QString, QMetaType::Int,    2,    2,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 5,    2,

       0        // eod
};

void FaceShiftClient::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        FaceShiftClient *_t = static_cast<FaceShiftClient *>(_o);
        switch (_id) {
        case 0: _t->askForConnection((*reinterpret_cast< QString(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 1: _t->initializeConnection(); break;
        case 2: _t->receivedNetworkUpdates((*reinterpret_cast< QSharedPointer<fs::fsMsgTrackingState>(*)>(_a[1]))); break;
        default: ;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (FaceShiftClient::*_t)(QString , int );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&FaceShiftClient::askForConnection)) {
                *result = 0;
            }
        }
    }
}

const QMetaObject FaceShiftClient::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_FaceShiftClient.data,
      qt_meta_data_FaceShiftClient,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *FaceShiftClient::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *FaceShiftClient::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_FaceShiftClient.stringdata))
        return static_cast<void*>(const_cast< FaceShiftClient*>(this));
    return QObject::qt_metacast(_clname);
}

int FaceShiftClient::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 3)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 3;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 3)
            *reinterpret_cast<int*>(_a[0]) = -1;
        _id -= 3;
    }
    return _id;
}

// SIGNAL 0
void FaceShiftClient::askForConnection(QString _t1, int _t2)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
struct qt_meta_stringdata_fs__FSTCPClient_t {
    QByteArrayData data[16];
    char stringdata[230];
};
#define QT_MOC_LITERAL(idx, ofs, len) \
    Q_STATIC_BYTE_ARRAY_DATA_HEADER_INITIALIZER_WITH_OFFSET(len, \
    qptrdiff(offsetof(qt_meta_stringdata_fs__FSTCPClient_t, stringdata) + ofs \
        - idx * sizeof(QByteArrayData)) \
    )
static const qt_meta_stringdata_fs__FSTCPClient_t qt_meta_stringdata_fs__FSTCPClient = {
    {
QT_MOC_LITERAL(0, 0, 15), // "fs::FSTCPClient"
QT_MOC_LITERAL(1, 16, 22), // "receivedNetworkUpdates"
QT_MOC_LITERAL(2, 39, 0), // ""
QT_MOC_LITERAL(3, 40, 38), // "QSharedPointer<fs::fsMsgTrack..."
QT_MOC_LITERAL(4, 79, 4), // "init"
QT_MOC_LITERAL(5, 84, 15), // "connectToServer"
QT_MOC_LITERAL(6, 100, 8), // "hostname"
QT_MOC_LITERAL(7, 109, 4), // "port"
QT_MOC_LITERAL(8, 114, 9), // "readyRead"
QT_MOC_LITERAL(9, 124, 5), // "error"
QT_MOC_LITERAL(10, 130, 28), // "QAbstractSocket::SocketError"
QT_MOC_LITERAL(11, 159, 12), // "stateChanged"
QT_MOC_LITERAL(12, 172, 28), // "QAbstractSocket::SocketState"
QT_MOC_LITERAL(13, 201, 5), // "state"
QT_MOC_LITERAL(14, 207, 9), // "connected"
QT_MOC_LITERAL(15, 217, 12) // "disconnected"

    },
    "fs::FSTCPClient\0receivedNetworkUpdates\0"
    "\0QSharedPointer<fs::fsMsgTrackingState>\0"
    "init\0connectToServer\0hostname\0port\0"
    "readyRead\0error\0QAbstractSocket::SocketError\0"
    "stateChanged\0QAbstractSocket::SocketState\0"
    "state\0connected\0disconnected"
};
#undef QT_MOC_LITERAL

static const uint qt_meta_data_fs__FSTCPClient[] = {

 // content:
       7,       // revision
       0,       // classname
       0,    0, // classinfo
       8,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       1,       // signalCount

 // signals: name, argc, parameters, tag, flags
       1,    1,   54,    2, 0x06 /* Public */,

 // slots: name, argc, parameters, tag, flags
       4,    0,   57,    2, 0x0a /* Public */,
       5,    2,   58,    2, 0x0a /* Public */,
       8,    0,   63,    2, 0x0a /* Public */,
       9,    1,   64,    2, 0x0a /* Public */,
      11,    1,   67,    2, 0x0a /* Public */,
      14,    0,   70,    2, 0x0a /* Public */,
      15,    0,   71,    2, 0x0a /* Public */,

 // signals: parameters
    QMetaType::Void, 0x80000000 | 3,    2,

 // slots: parameters
    QMetaType::Void,
    QMetaType::Void, QMetaType::QString, QMetaType::Int,    6,    7,
    QMetaType::Void,
    QMetaType::Void, 0x80000000 | 10,    2,
    QMetaType::Void, 0x80000000 | 12,   13,
    QMetaType::Void,
    QMetaType::Void,

       0        // eod
};

void fs::FSTCPClient::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        FSTCPClient *_t = static_cast<FSTCPClient *>(_o);
        switch (_id) {
        case 0: _t->receivedNetworkUpdates((*reinterpret_cast< QSharedPointer<fs::fsMsgTrackingState>(*)>(_a[1]))); break;
        case 1: _t->init(); break;
        case 2: _t->connectToServer((*reinterpret_cast< QString(*)>(_a[1])),(*reinterpret_cast< int(*)>(_a[2]))); break;
        case 3: _t->readyRead(); break;
        case 4: _t->error((*reinterpret_cast< QAbstractSocket::SocketError(*)>(_a[1]))); break;
        case 5: _t->stateChanged((*reinterpret_cast< QAbstractSocket::SocketState(*)>(_a[1]))); break;
        case 6: _t->connected(); break;
        case 7: _t->disconnected(); break;
        default: ;
        }
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        switch (_id) {
        default: *reinterpret_cast<int*>(_a[0]) = -1; break;
        case 4:
            switch (*reinterpret_cast<int*>(_a[1])) {
            default: *reinterpret_cast<int*>(_a[0]) = -1; break;
            case 0:
                *reinterpret_cast<int*>(_a[0]) = qRegisterMetaType< QAbstractSocket::SocketError >(); break;
            }
            break;
        case 5:
            switch (*reinterpret_cast<int*>(_a[1])) {
            default: *reinterpret_cast<int*>(_a[0]) = -1; break;
            case 0:
                *reinterpret_cast<int*>(_a[0]) = qRegisterMetaType< QAbstractSocket::SocketState >(); break;
            }
            break;
        }
    } else if (_c == QMetaObject::IndexOfMethod) {
        int *result = reinterpret_cast<int *>(_a[0]);
        void **func = reinterpret_cast<void **>(_a[1]);
        {
            typedef void (FSTCPClient::*_t)(QSharedPointer<fs::fsMsgTrackingState> );
            if (*reinterpret_cast<_t *>(func) == static_cast<_t>(&FSTCPClient::receivedNetworkUpdates)) {
                *result = 0;
            }
        }
    }
}

const QMetaObject fs::FSTCPClient::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_fs__FSTCPClient.data,
      qt_meta_data_fs__FSTCPClient,  qt_static_metacall, Q_NULLPTR, Q_NULLPTR}
};


const QMetaObject *fs::FSTCPClient::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->dynamicMetaObject() : &staticMetaObject;
}

void *fs::FSTCPClient::qt_metacast(const char *_clname)
{
    if (!_clname) return Q_NULLPTR;
    if (!strcmp(_clname, qt_meta_stringdata_fs__FSTCPClient.stringdata))
        return static_cast<void*>(const_cast< FSTCPClient*>(this));
    return QObject::qt_metacast(_clname);
}

int fs::FSTCPClient::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 8)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 8;
    } else if (_c == QMetaObject::RegisterMethodArgumentMetaType) {
        if (_id < 8)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 8;
    }
    return _id;
}

// SIGNAL 0
void fs::FSTCPClient::receivedNetworkUpdates(QSharedPointer<fs::fsMsgTrackingState> _t1)
{
    void *_a[] = { Q_NULLPTR, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 0, _a);
}
QT_END_MOC_NAMESPACE
