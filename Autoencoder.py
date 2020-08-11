from keras.models import Model  # 泛型模型
from keras.layers import Dense, Input
import pandas as pd

def encoder32(trainX):
    # 压缩特征维度至32维
    encoding_dim = 32

    # this is our input placeholder
    input_img = Input(shape=(70,))

    # 编码层
    encoded = Dense(70,activation='relu')(input_img)
    encoded = Dense(64, activation='relu')(encoded)
    encoder_output = Dense(encoding_dim)(encoded)

    # 解码层
    decoded = Dense(64, activation='relu')(encoder_output)
    decoded = Dense(70, activation='relu')(decoded)


    # 构建自编码模型
    autoencoder = Model(inputs=input_img, outputs=decoded)

    # 构建编码模型
    encoder = Model(inputs=input_img, outputs=encoder_output)

    # compile autoencoder
    autoencoder.compile(optimizer='adam', loss='mse')

    # training
    autoencoder.fit(trainX, trainX, epochs=6, batch_size=8, shuffle=True)

    # plotting
    encoded_imgs = encoder.predict(trainX)
    return encoded_imgs


def encoder64(trainX):
    # 压缩特征维度至64维
    encoding_dim = 64

    # this is our input placeholder
    input_img = Input(shape=(70,))

    # 编码层
    encoded = Dense(70, activation='relu')(input_img)
    encoder_output = Dense(encoding_dim)(encoded)

    # 解码层
    decoded = Dense(70, activation='relu')(encoder_output)

    # 构建自编码模型
    autoencoder = Model(inputs=input_img, outputs=decoded)

    # 构建编码模型
    encoder = Model(inputs=input_img, outputs=encoder_output)

    # compile autoencoder
    autoencoder.compile(optimizer='adam', loss='mse')

    # training
    autoencoder.fit(trainX, trainX, epochs=6, batch_size=8, shuffle=True)

    # plotting
    encoded_imgs = encoder.predict(trainX)
    return encoded_imgs

def getTemFeature(df):
    '''得到高级特征'''
    df = df.drop(0, axis=1)
    print(df.values.shape)

    train = df.values
    train = train.astype('float64')
    encoder3 = encoder32(train)
    encoder6 = encoder64(train)
    # print(encoder3.shape)
    # print(encoder6.shape)

    df_encoder32 = pd.DataFrame(encoder3)
    df_encoder64 = pd.DataFrame(encoder6)

    encoder = pd.concat([df_encoder32, df_encoder64], axis=1, ignore_index=True)
    return encoder